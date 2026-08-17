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
#   julia_csi_watref.raw same layout, water reference only
# Written to tmp_dir:
#   julia_recoinfo.m     CSI dimensions, read by julia_write_lcm_files on the MATLAB side

import Pkg
# MRSI.jl and JSON.jl may live in their own environment, see JULIA_MRSI_PKG in
# InstallProgramPaths.sh. If the variable is empty the default environment is used.
if get(ENV, "JULIA_MRSI_PKG", "") != ""
    Pkg.activate(ENV["JULIA_MRSI_PKG"]; io=devnull)
end
using JSON
using MRSI

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

if is_water_ref
    # "Method,/path/to/file" or just "Method", in which case the metabolite file is used
    watref_setting = get(p, "WaterReference_MethodAndFile", "")
    parts = split(watref_setting, ',')
    dat_file = length(parts) >= 2 ? strip(parts[2]) : csi_paths[1]
    isempty(dat_file) && (dat_file = csi_paths[1])
else
    dat_file = csi_paths[clamp(cur_avg, 1, n_files)]
end
isfile(dat_file) || error("CSI file not found: $dat_file")
println("Julia: reconstructing $dat_file (average $cur_avg, water reference: $is_water_ref)")

# Map the flags of InitialParameters to the keyword arguments of MRSI.reconstruct
hamming_flag = get(p, "hamming_flag", 0) == 1
noisedecor_flag = get(p, "noisedecorrelation_flag", 0) == 1
zero_fill_flag = get(p, "ZeroFillMetMaps_flag", 0) == 1

lipid_decon = nothing
if get(p, "LipidDecon_flag", 0) == 1
    method = uppercase(split(get(p, "LipidDecon_MethodAndNoOfLoops", "L2,10"), ',')[1])
    lipid_decon = method == "L1" ? :L1 : :L2
end

gradient_delay_us = [12.42 + 10.27im, 12.38 + 10.75im, 10.14 + 8.99im]  # MRSI.jl defaults
if get(p, "GradientDelay_flag", 0) == 1
    gd_str = string(get(p, "GradientDelay", ""))
    if !isempty(gd_str) && gd_str != "0"
        try
            gradient_delay_us = parse_gradient_delays(gd_str)
        catch e
            @warn "Could not read GradientDelay '$gd_str', using the MRSI.jl defaults: $e"
        end
    end
end

mmap_val = mmap_arg == "true" ? true : mmap_arg == "false" ? false : :auto

println("Julia: hamming=$hamming_flag, noise_decorrelation=$noisedecor_flag, lipid_decon=$lipid_decon")
println("Julia: gradient_delay_us=$gradient_delay_us, mmap=$mmap_val, zero_fill=$zero_fill_flag")

result = MRSI.reconstruct(
    dat_file;
    datatype = ComplexF32,
    mmap = mmap_val,
    do_noise_decorrelation = noisedecor_flag,
    do_hamming_filter = hamming_flag,
    do_hamming_filter_z = hamming_flag,
    lipid_decon = lipid_decon,
    gradient_delay_us = gradient_delay_us,
    zero_fill = zero_fill_flag,
)

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

if !is_water_ref && (n_files <= 1 || cur_avg == n_files)
    Nx, Ny = sz[1], sz[2]
    Nz = length(sz) >= 4 ? sz[3] : 1
    Nt = sz[end]
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
