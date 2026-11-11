#!/usr/bin/env julia
# Julia MRSI reconstruction entry script.
# Called by run_matlab.sh when the -S flag of Part1_ProcessMRSI.sh is used.
#
# Usage: julia run_julia_reco.jl <tmp_dir> <cur_avg> [mmap]
#
#   tmp_dir : temporary directory of the current run, holds InitialParameters.json
#   cur_avg : 1 based index of the current average
#   mmap    : "true", "false" or "auto" (memory mapping mode of MRSI.jl)
#
# The parameters are read from InitialParameters.json, which write_InitialParameters.sh
# writes next to InitialParameters.m. Do not parse the .m file here, it is MATLAB source.
#
# Written to out_path:
#   julia_csi.raw        accumulated complex float32, interleaved Re/Im, [Nx,Ny,Nz,Nt]
#   julia_csi_watref.raw same layout, water reference only, on its own grid
#   CombinedCSI.mat      csi.Data, the array MRSI_Reconstruction stores under the same
#                        name, read by julia_write_lcm_files on the MATLAB side
#   WaterReference.mat   weights.Data and WaterReferenceCSI, written on the water pass
#                        in the layout MRSI_Reconstruction.m uses, so either pipeline
#                        can read the other's coil weights; WaterReferenceCSI.Data is on
#                        the metabolite grid, and channel_mixing_matrix is the channel
#                        mixing both passes reconstruct with. A water reference with
#                        another slab (water_reference_grid.jl) leaves only the mixing
#                        on the water pass, and WaterReferenceUncombined.mat beside it;
#                        the metabolite pass combines that with its own coil weights,
#                        writes weights and WaterReferenceCSI, and removes it. No
#                        julia_csi_watref.raw then.
# Written to tmp_dir:
#   julia_recoinfo.m     CSI dimensions, read by julia_write_lcm_files on the MATLAB side

import Pkg
# MRSI.jl and JSON.jl may live in their own environment, see JULIA_MRSI_PKG in
# InstallProgramPaths.sh. If the variable is empty the default environment is used.
if get(ENV, "JULIA_MRSI_PKG", "") != ""
    Pkg.activate(ENV["JULIA_MRSI_PKG"]; io=devnull)
end
using JSON
using MAT
using MRSI
include(joinpath(@__DIR__, "water_reference_grid.jl"))

"""Parse a gradient delay string like "[12.56, 12.54, 10.08]" into complex microseconds.
The complex form "[12.42+10.27im, ...]" is also accepted."""
function parse_gradient_delays(s::AbstractString)
    s = strip(s, ['[', ']', ' '])
    result = ComplexF64[]
    for part in split(s, ',')
        part = strip(part)
        m = match(r"^(-?[\d.]+(?:[eE][+-]?\d+)?)\s*([+-]\s*[\d.]+(?:[eE][+-]?\d+)?)im$", part)
        if m !== nothing
            re = parse(Float64, replace(String(m[1]), " " => ""))
            im = parse(Float64, replace(String(m[2]), " " => ""))
            push!(result, complex(re, im))
        else
            push!(result, complex(parse(Float64, part), 0.0))
        end
    end
    return result
end

"""Write complex data as float32, interleaved Re/Im."""
function write_complex_raw(path::AbstractString, data::AbstractArray{<:Complex})
    open(path, "w") do fid
        for x in data
            write(fid, Float32(real(x)))
            write(fid, Float32(imag(x)))
        end
    end
end

# MRSI_Reconstruction.m reduces the reference to FID point 4 before taking the conj
# (`:566-570`, then `:958`), so weights written here use the same point to stay
# interchangeable with the ones it stores.
const MATLAB_REF_POINT = 4

"""How MRSI_Reconstruction.m combines a scan: MUSICAL at `:963`, or the phase-only
"Phase VC" branch at `:968` for a single channel. Deriving with the wrong one gives
weights that are simply different, not slightly off."""
combine_method_of(file) = MRSI.read_scan_info(file)[:ONLINE][:n_channels] == 1 ? :phase_only : :musical

"""WaterReference.mat in the layout MRSI_Reconstruction.m writes. `WaterReferenceCSI`
carries the coil-compression count the metabolite pass compares against its own, and the
water signal as W1 reads it: on the metabolite grid, scaled by `gain_ratio`. `dwelltime`
(ns) is that of `water_csi`, `par` any further entries of its Par."""
function write_water_reference(path, weights, water_csi, mixing, water_info, gain_ratio;
                               dwelltime=water_info[:ONLINE][:dwelltime], par=Dict{String,Any}())
    println("Julia: writing $path")
    # RecoPar.Dwelltimes lets julia_write_lcm_files refuse a water reference LCModel
    # would read with the wrong bandwidth.
    contents = Dict{String,Any}(
        "weights" => Dict("Data" => weights),
        "WaterReferenceCSI" => Dict(
            "Data" => water_csi,
            "Par" => merge(Dict{String,Any}("coilcompression_to_numb_coils" => water_info[:n_compressed_coils],
                                            "AcquiredMatrix" => float.(collect(matrix_size(water_info))),
                                            "ReconstructionGainRatio" => gain_ratio), par),
            "RecoPar" => Dict("Dwelltimes" => float(dwelltime),
                              "vecSize" => float(size(water_csi, ndims(water_csi))))))
    isnothing(mixing) || (contents["channel_mixing_matrix"] = mixing)
    matwrite(path, contents)
end

"""The channel mixing a water pass stored with its weights, or `nothing`."""
function stored_channel_mixing(path)
    isfile(path) || return nothing
    mixing = get(matread(path), "channel_mixing_matrix", nothing)
    isnothing(mixing) && return nothing
    # matread returns a 1x1 matrix as a scalar
    return mixing isa Number ? fill(ComplexF64(mixing), 1, 1) : Matrix{ComplexF64}(mixing)
end

"""The water file of `-w`, the last comma-separated field, as GetPar_CreateTempl_MaskPart1.m
reads it: "W1,/path/to/file", or a bare path whose name carries the method."""
water_reference_file(p) = String(strip(last(split(get(p, "WaterReference_MethodAndFile", ""), ','))))

"""(nx, ny, nz) of the reconstructed CSI."""
matrix_size(scan_info) =
    (scan_info[:ONLINE][:n_frequency], scan_info[:ONLINE][:n_frequency], scan_info[:ONLINE][:n_part])

# How far two FOVs may differ for a k-space resize between them
const FOV_TOLERANCE_MM = 0.1
const FOV_TOLERANCE_RAD = 1e-3

"""How the FOVs of two scans differ, as text; empty when they are the same."""
function fov_differences(met_info, water_info)
    differences = String[]
    for key in (:fov_readout, :fov_phase, :fov_slice)
        abs(met_info[key] - water_info[key]) <= FOV_TOLERANCE_MM ||
            push!(differences, "$key $(met_info[key]) vs $(water_info[key]) mm")
    end
    shift = sqrt(sum(abs2, met_info[:position] .- water_info[:position]))
    shift <= FOV_TOLERANCE_MM || push!(differences, "centres $(round(shift; digits=2)) mm apart")
    n1, n2 = met_info[:slice_normal], water_info[:slice_normal]
    tilt = acos(clamp(sum(n1 .* n2) / sqrt(sum(abs2, n1) * sum(abs2, n2)), -1, 1))
    tilt <= FOV_TOLERANCE_RAD || push!(differences, "slab normals $(round(rad2deg(tilt); digits=3)) deg apart")
    rotation = abs(rem2pi(met_info[:in_plane_rotation] - water_info[:in_plane_rotation], RoundNearest))
    rotation <= FOV_TOLERANCE_RAD ||
        push!(differences, "in-plane rotations $(round(rad2deg(rotation); digits=3)) deg apart")
    return differences
end

"""Refuse a water reference that does not cover the metabolite scan's FOV. Zero-filling
or cutting k-space moves data between matrices of one FOV only; another slab takes the
route of water_reference_grid.jl instead."""
function check_same_fov(met_info, water_info, met_file, water_file)
    differences = fov_differences(met_info, water_info)
    isempty(differences) && return
    error("The water reference $water_file does not share the FOV of $met_file: " *
          join(differences, ", ") * ". Its coil weights and water signal can only be moved " *
          "between matrices of the same FOV.")
end

"""How the water pass brings the water onto the metabolite grid: `:resize` in k-space
when both scans share the FOV, `:metabolite_grid` when only the slab differs. Any other
difference stays with `:resize`, which refuses it."""
function water_route(met_info, water_info)
    isempty(fov_differences(met_info, water_info)) && return :resize
    slab_only = try
        slab_offset_mm(met_info, water_info)
        true
    catch
        false
    end
    return slab_only ? :metabolite_grid : :resize
end

"""`data` on the water reference's grid, brought onto the metabolite scan's by
zero-filling or cutting k-space, as MRSI_Reconstruction.m's op_CoilCombineData does
with the weights it loads."""
function onto_metabolite_grid(data, met_info, water_info, met_file, water_file)
    grid = matrix_size(met_info)
    size(data)[1:3] == grid && return data
    check_same_fov(met_info, water_info, met_file, water_file)
    println("Julia: resizing $(size(data)[1:3]) to the metabolite grid $grid in k-space")
    return MRSI.resize_grid(data, grid)
end

"""Metabolite over water reconstruction gain. With both scans combined with the same
coil weights neither is normalised by its own reference scan, so the images of one
density differ by this ratio, and W1 has to scale the water signal by it."""
reconstruction_gain_ratio(met_info, water_info; kw...) =
    MRSI.reconstruction_gain(met_info[:ONLINE]; kw...) / MRSI.reconstruction_gain(water_info[:ONLINE]; kw...)

"""The water pass for a water reference with another slab: the water uncombined on the
metabolite grid (water_reference_grid.jl), mixed as the metabolite scan, in `path`. Then
WaterReference.mat with the mixing alone, last, since its existence marks the pass done."""
function write_uncombined_water(path, watref_path, water_scan, met_scan, mixing; do_hamming_filter, gradient_delay_us)
    water, met = water_scan[:ONLINE], met_scan[:ONLINE]
    offset = slab_offset_mm(met, water)
    println("Julia: water slab $(water[:fov_slice]) mm in $(water[:n_part]) partitions, metabolite slab " *
            "$(met[:fov_slice]) mm in $(met[:n_part]), centres $(round(offset; digits=3)) mm apart along the normal")
    data = water_onto_metabolite_grid(water_scan, met_scan; n_points=water[:n_fid], channel_mixing_matrix=mixing,
                                      do_hamming_filter, gradient_delay_us)
    println("Julia: writing $path")
    matwrite(path, Dict{String,Any}("Data" => data, "SlabOffset_mm" => offset))
    println("Julia: writing $watref_path")
    matwrite(watref_path, Dict{String,Any}("channel_mixing_matrix" => mixing))
end

"""The metabolite pass's half: the water left by `write_uncombined_water`, combined with
the coil weights the metabolite CSI was combined with (ICE: its PATREFSCAN `patref` at
point 5), so the receive field cancels and no gain is left between the two; resampled to
the metabolite dwell time, at most as long as the metabolite FID, and written as
WaterReferenceCSI with those weights and the metabolite scan's channel count."""
function combine_uncombined_water(path, watref_path, patref, mixing, water_scan, met_scan)
    stored = matread(path)
    weights = MRSI.coil_combine_weights(patref)
    water = MRSI.apply_coil_weights_and_combine(stored["Data"], weights)
    water_dwell, met_dwell = water_scan[:ONLINE][:dwelltime], met_scan[:ONLINE][:dwelltime]
    println("Julia: water reference combined with the metabolite scan's coil weights, resampled from " *
            "$(water_dwell / 1e3) to $(met_dwell / 1e3) us")
    water = resample_fid(water, water_dwell, met_dwell)
    n_met = met_scan[:ONLINE][:n_fid]
    size(water, 4) > n_met && (water = water[:, :, :, 1:n_met])
    par = Dict{String,Any}("Route" => "metabolite_grid", "SlabOffset_mm" => stored["SlabOffset_mm"],
                           "AcquiredDwelltime" => float(water_dwell), "CoilWeights" => "metabolite PATREFSCAN",
                           "coilcompression_to_numb_coils" => met_scan[:n_compressed_coils])
    write_water_reference(watref_path, weights, water, mixing, water_scan, 1.0; dwelltime=met_dwell, par)
    rm(path)
end

"""What deep_crt_mrsi.combined_csi reads out of a CombinedCSI.mat.

It wants three things: the combined spectra as `csi.Data`, the UNCOMBINED water
reference as `image_FullFID.Data`, and the acquisition parameters under
`csi.RecoPar`. Writing only the spectra leaves the deep fitting unable to run at
all, which is why the field names here are its names and not ours."""
function combined_csi_contents(csi_data, dat_file, recon_kw, patref_override=nothing)
    contents = Dict{String,Any}("csi" => Dict{String,Any}("Data" => csi_data))
    scan_info = MRSI.read_scan_info(dat_file)
    online = scan_info[:ONLINE]

    contents["csi"]["RecoPar"] = Dict{String,Any}(
        "Dwelltimes" => float(online[:dwelltime]),            # ns, as combined_csi expects
        "LarmorFreq" => float(online[:larmor_frequency]),
        "FoV_Read" => float(online[:fov_readout]),
        "FoV_Partition" => float(online[:fov_slice]),
        "nFreqEnc_FinalMatrix" => float(size(csi_data, 1)),   # after reconstruction, not the encoded size
        # water_scaling.py's tissue model takes the PSF from it
        "KSpaceWindow" => String(get(recon_kw, :kspace_window, MRSI.default_kspace_window(scan_info))),
    )

    # MindGrab (-m mindgrab) needs the affine of the CSI grid. An MRSI.jl without
    # slice_aligned_rotation pairs the in-plane axes the wrong way round.
    if isdefined(MRSI, :slice_aligned_rotation)
        contents["csi"]["RecoPar"]["PatientAffineLPS"] =
            Matrix{Float64}(MRSI.build_patient_affine(scan_info; axis_convention=:lps))
    end

    patref = isnothing(patref_override) ? uncombined_patref(scan_info, recon_kw) :
                                          as_uncombined_5d(patref_override)
    if isnothing(patref)
        @warn "No PATREFSCAN, so CombinedCSI.mat carries no image_FullFID and the deep fitting cannot read it"
    else
        contents["image_FullFID"] = Dict{String,Any}("Data" => patref)
    end
    return contents
end

"""Whatever the reconstruction handed back, as (x, y, z, t, cha)."""
function as_uncombined_5d(patref)
    patref = patref isa AbstractVector ? first(patref) : patref
    data = Array{ComplexF32}(patref)
    return ndims(data) == 5 ? data :
           reshape(data, size(data, 1), size(data, 2), size(data, 3), size(data, 4), :)
end

"""The water reference with its channels intact, as (x, y, z, t, cha).

The deep fitting does its own coil combination, so it needs the uncombined
reference; a combined one would be silently accepted as single-channel."""
function uncombined_patref(scan_info, recon_kw)
    haskey(scan_info, :PATREFSCAN) || return nothing
    reco = MRSI.reconstruct_uncombined(scan_info; scans=[:PATREFSCAN], recon_kw...)
    patref = get(reco, :PATREFSCAN, nothing)
    isnothing(patref) && return nothing
    patref = patref isa AbstractVector ? first(patref) : patref
    data = Array{ComplexF32}(patref)
    return ndims(data) == 5 ? data : reshape(data, size(data, 1), size(data, 2), size(data, 3), size(data, 4), :)
end

"""`weights.Data` from a WaterReference.mat, whichever pipeline wrote it, or `nothing`.

MATLAB drops trailing singleton dimensions when it saves, so a single-channel file
comes back as (x, y, z) rather than the (x, y, z, 1, cha) the combination indexes."""
function stored_coil_weights(path)
    isfile(path) || return nothing
    stored = matread(path)
    weights = get(stored, "weights", nothing)
    weights isa AbstractDict && haskey(weights, "Data") || return nothing
    w = weights["Data"]
    return ndims(w) == 5 ? w : reshape(w, size(w, 1), size(w, 2), size(w, 3), 1, :)
end

"""Read a file written by write_complex_raw back into an array of size sz."""
function read_complex_raw(path::AbstractString, sz::NTuple)
    n = prod(sz)
    buf = Array{Float32}(undef, 2 * n)
    read!(path, buf)
    data = Array{ComplexF32}(undef, sz...)
    for i in 1:n
        data[i] = ComplexF32(buf[2i-1], buf[2i])
    end
    return data
end

# Arguments
length(ARGS) < 2 && error("Usage: julia run_julia_reco.jl <tmp_dir> <cur_avg> [mmap]")
tmp_dir = ARGS[1]
cur_avg = parse(Int, ARGS[2])
mmap_arg = length(ARGS) >= 3 ? ARGS[3] : "auto"
# Which reconstruction this one reproduces. "ice" is the online route the scanner
# runs; "matlab" is MRSI_Reconstruction.m. The two disagree about the frequency
# drift, the trajectory delays and several smaller choices, so a run has to say
# which of them it is trying to be rather than land between the two.
parity_mode = length(ARGS) >= 4 && !isempty(ARGS[4]) ? lowercase(ARGS[4]) : "ice"
parity_mode in ("ice", "matlab") ||
    error("Unknown -S mode '$(parity_mode)'. Use ice or matlab.")

par_file = joinpath(tmp_dir, "InitialParameters.json")
isfile(par_file) || error("InitialParameters.json not found in $tmp_dir")
p = JSON.parsefile(par_file)

out_path = get(p, "out_path", "")
csi_paths = get(p, "csi_path", String[])
n_files = length(csi_paths)
isempty(out_path) && error("out_path missing in $par_file")
isempty(csi_paths) && error("csi_path missing in $par_file")

# The water reference is reconstructed first, before WaterReference.mat exists
is_water_ref = false
if get(p, "WaterReference_flag", 0) == 1 && !isfile(joinpath(out_path, "WaterReference.mat"))
    is_water_ref = true
    println("Julia: water reference pass.")
end

# GetPar has already checked that the water file exists.
water_file = get(p, "WaterReference_flag", 0) == 1 ? water_reference_file(p) : nothing
met_file = csi_paths[1]
dat_file = is_water_ref ? water_file : csi_paths[clamp(cur_avg, 1, n_files)]
isfile(dat_file) || error("CSI file not found: $dat_file")
println("Julia: reconstructing $dat_file (average $cur_avg, water reference: $is_water_ref)")

# Map the flags of InitialParameters to the keyword arguments of MRSI.reconstruct
hamming_flag = get(p, "hamming_flag", 0) == 1
noisedecor_flag = get(p, "noisedecorrelation_flag", 0) == 1
zero_fill_flag = get(p, "ZeroFillMetMaps_flag", 0) == 1

lipid_decon = nothing
lipid_decon_kw = NamedTuple()
if get(p, "LipidDecon_flag", 0) == 1
    spec = split(get(p, "LipidDecon_MethodAndNoOfLoops", "L2,10"), ',')
    method = uppercase(spec[1])
    value = length(spec) > 1 ? tryparse(Float64, spec[2]) : nothing
    # WALRUS (formerly WALINET) is the alternative to the regularization, not an
    # addition to it. It runs after the reconstruction, in walrus_clean_csi.py, so
    # nothing happens here.
    lipid_decon = method in ("WALRUS", "WALINET") ? nothing : method == "L1" ? :L1 : :L2

    if lipid_decon == :L1 && !isnothing(value)
        # An iteration count means the same thing on both sides.
        lipid_decon_kw = (; L1_n_loops = round(Int, value))
    elseif lipid_decon == :L2 && !isnothing(value)
        # It does not follow for L2. MATLAB takes LipidDecon_L2BetaCorrFactor, a
        # multiplier on a beta it computes, default 1; MRSI.jl takes L2_beta, an
        # absolute strength, default 0.2. The same number means different things,
        # so it is not forwarded, and saying so beats discarding it in silence.
        @warn "-L \"$(method),$(spec[2])\": the Julia reconstruction uses MRSI.jl's " *
              "own L2 strength (L2_beta=0.2). The value is MATLAB's correction " *
              "factor and is not on that scale, so it is not applied here."
    end
end

# -G is the only thing that asks for a trajectory delay here. Part1 documents
# "If -G is not used no gradient delay is used", which is what the MATLAB
# reconstruction does; the online route instead takes the delays whenever the
# protocol asks for them. So the value is only pinned when the caller gave one,
# and each mode fills the gap the way its reference does.
gradient_delay_us = nothing
if get(p, "GradientDelay_flag", 0) == 1
    gd_str = string(get(p, "GradientDelay", ""))
    if !isempty(gd_str) && gd_str != "0"
        try
            # "global" is required: try opens a soft scope, so a bare assignment
            # would create a local and the parsed delays would be discarded.
            global gradient_delay_us = parse_gradient_delays(gd_str)
        catch e
            @warn "Could not read GradientDelay '$gd_str', reconstructing without a delay: $e"
        end
    end
end

mmap_val = mmap_arg == "true" ? true : mmap_arg == "false" ? false : :auto

println("Julia: hamming=$hamming_flag, noise_decorrelation=$noisedecor_flag, lipid_decon=$lipid_decon")
println("Julia: mode=$parity_mode, gradient_delay_us=$(something(gradient_delay_us, "protocol")), " *
        "mmap=$mmap_val, zero_fill=$zero_fill_flag")

# Coil weights fit only the channels they were derived in, and coil compression and
# noise decorrelation are derived per measurement. So both passes mix the channels as
# the metabolite scan does, and the water pass stores that mixing with its weights.
watref_path = joinpath(out_path, "WaterReference.mat")
shared_mixing = if is_water_ref
    Matrix{ComplexF64}(MRSI.channel_mixing_matrix(met_file; do_noise_decorrelation=noisedecor_flag))
else
    stored_channel_mixing(watref_path)
end
isnothing(shared_mixing) ||
    println("Julia: channels mixed as in $met_file, $(size(shared_mixing, 2)) channels")
mixing_kw = isnothing(shared_mixing) ? (;) : (; channel_mixing_matrix = shared_mixing)

# What Part1's own flags ask for, in both modes.
common_kw = (
    datatype = ComplexF32,
    mmap = mmap_val,
    do_noise_decorrelation = noisedecor_flag,
    do_hamming_filter = hamming_flag,
    lipid_decon = lipid_decon,
    zero_fill = zero_fill_flag,
    lipid_decon_kw...,
    mixing_kw...,
)

# The settings MATLAB picks differently. Each one was measured against a
# MRSI_Reconstruction.m run of the same data: with the delays and the drift fit
# matched, the two reference scans agree exactly.
matlab_kw = (
    ice = false,
    do_hamming_filter_z = hamming_flag ? :matlab : false,
    frequency_drift_correction = nothing,
    ref_point_for_combine = MATLAB_REF_POINT,
    gradient_delay_us = something(gradient_delay_us, [0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im]),
    # op_ReconstructNonCartMRData.m's rotation of the FOV offset into the slice
    fov_shift_rotation = :matlab_legacy,
    # MATLAB windows with Hamming whatever the protocol states
    kspace_window = :hamming,
)

# -G, when it was given, applies in both modes. Left out otherwise so the ICE
# mode's protocol gating decides.
delay_kw = isnothing(gradient_delay_us) ? (;) : (; gradient_delay_us)

# The ICE mode shifts the FOV in the scanner's slice basis.
ice_kw = (; fov_shift_rotation = :siemens, delay_kw...)

# An MRSI.jl without the rotation modes ignores the keyword and shifts with its own.
isdefined(MRSI, :matlab_fov_shift_rotation) ||
    error("This MRSI.jl has no fov_shift_rotation=:matlab_legacy; JULIA_MRSI_PKG needs a newer one.")

# What the helpers reconstruct the reference scan and the coil weights with. It
# has to be the settings the CSI was reconstructed with, or the field map and the
# weights come from a differently reconstructed reference than the data.
recon_kw = parity_mode == "matlab" ? (; common_kw..., matlab_kw...) :
                                     (; common_kw..., ice_kw...)

# A water reference with another slab is reconstructed onto the metabolite grid here and
# combined in the metabolite pass, with that scan's own weights.
uncombined_path = joinpath(out_path, "WaterReferenceUncombined.mat")
# The water pass takes the metabolite scan's k-space window, so both share one PSF.
# The MATLAB mode has Hamming in matlab_kw.
water_window_kw = (;)
if is_water_ref
    met_scan, water_scan = MRSI.read_scan_info(met_file), MRSI.read_scan_info(dat_file)
    parity_mode == "ice" && (water_window_kw = (; kspace_window = MRSI.default_kspace_window(met_scan)))
    if water_route(met_scan[:ONLINE], water_scan[:ONLINE]) === :metabolite_grid
        parity_mode == "ice" ||
            error("The water reference $dat_file has another slab than $met_file; that route " *
                  "reconstructs in ICE mode only.")
        write_uncombined_water(uncombined_path, watref_path, water_scan, met_scan, shared_mixing;
                               do_hamming_filter=hamming_flag, gradient_delay_us)
        println("Julia: done.")
        exit(0)
    end
end

# The water reference supplies the coil weights for both scans, the way
# MRSI_Reconstruction.m reuses the ones it stored in WaterReference.mat. The water
# pass stores the weights its own CSI is combined with, at MATLAB's FID point, so
# the ratio of the two CSIs is the ratio of their signals, up to the reconstruction
# gains. The ICE route takes them from its reconstruction, further down.
coil_weights = nothing
if is_water_ref && parity_mode == "matlab"
    # Derived with the settings the CSI is reconstructed with, so they are the ones
    # its combination would derive itself.
    coil_weights = MRSI.coil_combine_weights(dat_file; ref_point_for_combine=MATLAB_REF_POINT,
                                             combine_method=combine_method_of(dat_file), recon_kw...)
elseif !is_water_ref
    coil_weights = stored_coil_weights(watref_path)
end
if !is_water_ref && !isnothing(coil_weights)
    met_info = MRSI.read_scan_info(dat_file)
    if size(coil_weights)[1:3] != matrix_size(met_info)
        isnothing(water_file) &&
            error("The coil weights in $watref_path are on another grid, and without -w there is " *
                  "no water reference to check its FOV against.")
        coil_weights = onto_metabolite_grid(coil_weights, met_info, MRSI.read_scan_info(water_file),
                                            dat_file, water_file)
    end
end
isnothing(coil_weights) || println("Julia: combining with the water reference coil weights")

ice_patref = nothing
result = try
    if parity_mode == "matlab"
        MRSI.reconstruct(dat_file; common_kw..., matlab_kw..., combine_method=combine_method_of(dat_file), coil_weights)
    else
        # reconstruct_like_ICE is the entry point ICE parity was established against:
        # it reads alICEProgramPara for the frequency offset, the drift and the
        # trajectory delays, exactly as ICE's Configurator does. Calling it is what
        # makes this mode ICE-equivalent; mirroring its choices here would be a
        # second copy to drift out of sync. prescan=false because Part1 has no use
        # for the reconstructed PATREFANDIMASCAN and it is not free.
        water_kw = is_water_ref ? (; ref_point_for_combine=MATLAB_REF_POINT,
                                     combine_method=combine_method_of(dat_file)) : (;)
        # combine=true: Part1 takes coil-combined spectra only, whatever the protocol asks.
        reco = MRSI.reconstruct_like_ICE(dat_file; prescan=false, common_kw..., ice_kw..., water_kw...,
                                         water_window_kw..., coil_weights, combine=true)
        # Keep the reference scan it already reconstructed: reconstruct_uncombined
        # would redo it through plain reconstruct, without the protocol gating, and
        # the field map would then come from a different reference than the data.
        global ice_patref = get(reco, :PATREFSCAN_UNCOMBINED, nothing)
        if is_water_ref
            isnothing(ice_patref) && error("The water reference $dat_file has no PATREFSCAN, so it cannot supply coil weights")
            global coil_weights = MRSI.coil_combine_weights(first(ice_patref); water_kw...)
        end
        reco[:ONLINE]
    end
catch err
    # The water reference's weights are per voxel of its own grid and per channel of
    # its own compression. A DimensionMismatch deep in the combination says nothing
    # about either, so the context goes first and the error keeps its trace.
    if !isnothing(coil_weights) && err isa DimensionMismatch
        @error "The coil weights from $watref_path have size $(size(coil_weights)) " *
               "(grid, 1, channels). Weights derived with another coil compression " *
               "cannot combine $dat_file."
    end
    rethrow()
end

# MRSI.reconstruct returns the CSI array of the single repetition, but the Vector of
# all repetitions as soon as the dataset has more than one.
csi_data = if result isa AbstractVector{<:AbstractArray}
    length(result) == 1 || error("$dat_file has $(length(result)) repetitions. " *
                                 "This pipeline averages over files, not over repetitions, " *
                                 "so the Julia reconstruction cannot be used for it.")
    Array{ComplexF32}(result[1])
elseif result isa AbstractArray
    Array{ComplexF32}(result)
else
    error("Unexpected return type of MRSI.reconstruct(): $(typeof(result))")
end

sz = size(csi_data)
println("Julia: reconstruction done, CSI size $sz")

out_raw = joinpath(out_path, is_water_ref ? "julia_csi_watref.raw" : "julia_csi.raw")

# Average the metabolite files the same way MRSI_Reconstruction.m does: add up
# over the calls and divide by the number of files on the last one.
if !is_water_ref && n_files > 1
    if cur_avg > 1 && isfile(out_raw)
        println("Julia: adding average $cur_avg to $out_raw")
        csi_data .+= read_complex_raw(out_raw, sz)
    end
    if cur_avg == n_files
        println("Julia: dividing by $n_files averages")
        csi_data ./= Float32(n_files)
    end
end

println("Julia: writing $out_raw")
write_complex_raw(out_raw, csi_data)

# W1 reads the water signal voxel by voxel against the spectra, so it is stored on
# their grid and in the metabolite reconstruction's gain.
if is_water_ref
    println("Julia: storing the coil weights of FID point $MATLAB_REF_POINT")
    met_info, water_info = MRSI.read_scan_info(met_file), MRSI.read_scan_info(dat_file)
    # reconstruct_like_ICE sets ice and its own z window
    gain_kw = parity_mode == "matlab" ? (; ice = false, matlab_kw.do_hamming_filter_z, matlab_kw.kspace_window) :
                                        (; ice = true, do_hamming_filter_z = :ice, water_window_kw...)
    gain_ratio = reconstruction_gain_ratio(met_info, water_info; gain_kw...)
    println("Julia: metabolite / water reconstruction gain $(round(gain_ratio; sigdigits=5)), " *
            "applied to the water signal")
    water_for_w1 = onto_metabolite_grid(csi_data, met_info, water_info, met_file, dat_file) .* Float32(gain_ratio)
    write_water_reference(watref_path, coil_weights, water_for_w1, shared_mixing, water_info, gain_ratio)
end

if !is_water_ref && isfile(uncombined_path)
    isnothing(ice_patref) &&
        error("$dat_file has no PATREFSCAN, whose coil weights combine the water reference on its grid.")
    combine_uncombined_water(uncombined_path, watref_path, first(ice_patref), shared_mixing,
                             MRSI.read_scan_info(water_file), MRSI.read_scan_info(dat_file))
end

if !is_water_ref && (n_files <= 1 || cur_avg == n_files)
    Nx, Ny = sz[1], sz[2]
    Nz = length(sz) >= 4 ? sz[3] : 1
    Nt = sz[end]
    combined_path = joinpath(out_path, "CombinedCSI.mat")
    println("Julia: writing $combined_path")
    matwrite(combined_path, combined_csi_contents(csi_data, dat_file, recon_kw, ice_patref))

    recoinfo_path = joinpath(tmp_dir, "julia_recoinfo.m")
    println("Julia: writing $recoinfo_path")
    open(recoinfo_path, "w") do f
        println(f, "% Written by run_julia_reco.jl, do not edit")
        println(f, "julia_csi_Nx = $Nx;")
        println(f, "julia_csi_Ny = $Ny;")
        println(f, "julia_csi_Nz = $Nz;")
        println(f, "julia_csi_Nt = $Nt;")
        println(f, "julia_csi_raw = '$(joinpath(out_path, "julia_csi.raw"))';")
        println(f, "julia_csi_watref_raw = '$(joinpath(out_path, "julia_csi_watref.raw"))';")
    end
end

println("Julia: done.")
