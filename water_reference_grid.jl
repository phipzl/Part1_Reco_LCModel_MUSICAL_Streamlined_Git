# A water reference acquired with another slab than the metabolite scan, brought onto the
# metabolite scan's grid with the metabolite reconstruction's point spread function.
# Included by run_julia_reco.jl and by tests/test_water_reference_grid.jl.
#
# In-plane both scans share the FOV and the CRT circles, so the water's circles up to the
# metabolite's outermost one are density compensated and windowed as the metabolite
# reconstruction does it and transformed onto its matrix. Along the slab each water
# partition is a sample of the object's kz spectrum; the metabolite partitions are
# evaluated from those samples inside the metabolite's kz extent, with its z window, at
# their positions in the water slab's frame. The result is the uncombined water as the
# metabolite reconstruction would have produced it, so combining it with the metabolite
# scan's own coil weights cancels the receive field and leaves no gain between the scans.

using MRSI
using LinearAlgebra

# Geometry tolerances, as run_julia_reco.jl's check_same_fov
const GRID_TOLERANCE_MM = 0.1
const GRID_TOLERANCE_RAD = 1e-3

"""Unit slab normal of a scan."""
slab_normal(info) = normalize(float.(info[:slice_normal]))

"""
    slab_offset_mm(met_info, water_info)

How far the water slab's centre lies from the metabolite slab's along their shared normal,
in mm, positive along the normal. An error unless the two scans share the in-plane FOV,
its centre, the normal and the in-plane rotation, which is what this route assumes.
"""
function slab_offset_mm(met_info, water_info)
    differences = String[]
    for key in (:fov_readout, :fov_phase)
        abs(met_info[key] - water_info[key]) <= GRID_TOLERANCE_MM ||
            push!(differences, "$key $(met_info[key]) vs $(water_info[key]) mm")
    end
    n_met, n_water = slab_normal(met_info), slab_normal(water_info)
    tilt = acos(clamp(dot(n_met, n_water), -1, 1))
    tilt <= GRID_TOLERANCE_RAD || push!(differences, "slab normals $(round(rad2deg(tilt); digits=3)) deg apart")
    rotation = abs(rem2pi(met_info[:in_plane_rotation] - water_info[:in_plane_rotation], RoundNearest))
    rotation <= GRID_TOLERANCE_RAD ||
        push!(differences, "in-plane rotations $(round(rad2deg(rotation); digits=3)) deg apart")
    shift = float.(water_info[:position]) .- float.(met_info[:position])
    along = dot(shift, n_met)
    in_plane = norm(shift .- along .* n_met)
    in_plane <= GRID_TOLERANCE_MM ||
        push!(differences, "FOV centres $(round(in_plane; digits=2)) mm apart in-plane")
    isempty(differences) ||
        error("The water reference's slab can only differ from the metabolite scan's along the normal: " *
              join(differences, ", ") * ".")
    return along
end

"""
    MetaboliteGrid(met_info, water_info)

What the water reconstruction takes from the metabolite scan: its matrix, the density
compensation of each circle, the circles of each partition, the z window and where its
partitions lie in the water slab's frame (mm from the water slab's centre along the normal,
the direction the partition index grows in).
"""
struct MetaboliteGrid
    n_grid::Int
    n_part::Int
    k_scale::Float64            # water's normalised k-space coordinates to the metabolite grid's
    area::Vector{Float64}       # effective area of each metabolite circle
    radii::Vector{Float64}      # radius of each metabolite circle as acquired, normalised
    circles::Vector{Int}        # outermost circle of each metabolite partition
    z_window::Bool
    fov_slice::Float64
    positions::Vector{Float64}  # metabolite partition centres, mm from the water slab centre
end

function MetaboliteGrid(met_info, water_info; offset_mm=slab_offset_mm(met_info, water_info))
    n_part = met_info[:n_part]
    circles = zeros(Int, n_part)
    for c in MRSI.sort_into_circles(met_info, 1)
        circles[c[:part]] = max(circles[c[:part]], c[:circle])
    end
    # the metabolite grid's partition q lies (q - (n÷2+1)) slices from its own centre, as
    # fft_slice_dim puts it and build_patient_affine describes it
    centre = n_part ÷ 2 + 1
    dz = met_info[:fov_slice] / n_part
    positions = ((1:n_part) .- centre) .* dz .- offset_mm
    return MetaboliteGrid(met_info[:n_frequency], n_part, MRSI.max_r(water_info) / MRSI.max_r(met_info),
                          MRSI.calculate_effective_area_per_circle(MRSI.get_radii(met_info; ice=true)),
                          collect(float.(met_info[:radii])), circles, MRSI.applies_z_apodization(:ice, n_part), met_info[:fov_slice],
                          collect(float.(positions)))
end

"""Metabolite partition index (from -n÷2) a kz in 1/mm falls into, and the continuous one."""
function metabolite_partition(grid::MetaboliteGrid, kz)
    p = kz * grid.fov_slice
    return round(Int, p), p
end

"""
Whether the metabolite scan samples the water circle `circle` at `kz`: the metabolite partition
cell kz falls into acquires a circle with that index.
"""
function in_metabolite_extent(grid::MetaboliteGrid, circle, kz)
    n, _ = metabolite_partition(grid, kz)
    half = grid.n_part ÷ 2
    abs(n) <= half || return false
    return circle <= grid.circles[n+half+1]
end

"""The metabolite's z window at kz, continued between its partitions."""
function z_weight(grid::MetaboliteGrid, kz)
    grid.z_window || return 1.0
    half = grid.n_part ÷ 2
    _, p = metabolite_partition(grid, kz)
    return MRSI.hamming_filter(clamp(p, -half, half) / (2half))
end

"""
    z_operator(grid, water_kz, fov_slice_water)

Matrix taking the water's per-partition in-plane images (kz = `water_kz`, 1/mm) to the
metabolite partitions: the partition FFT of the metabolite reconstruction, evaluated at the
metabolite partitions' positions, with its z window and its kz extent. A water partition
stands for 1/fov_slice_water of kz where a metabolite one stands for 1/fov_slice, hence the
ratio of the two slabs.
"""
function z_operator(grid::MetaboliteGrid, water_kz, fov_slice_water)
    E = zeros(ComplexF64, grid.n_part, length(water_kz))
    quadrature = grid.fov_slice / fov_slice_water
    for (j, kz) in enumerate(water_kz), q in 1:grid.n_part
        E[q, j] = quadrature * z_weight(grid, kz) * cispi(-2 * kz * grid.positions[q])
    end
    return E
end

"""
    water_kdata(c; channel_mixing_matrix, frequency_drift, gradient_delay_us, do_freq_cor, fov_shift_rotation)

One water circle's k-space data up to the density compensation, as `reconstruct(::Circle)`
prepares it: read and mixed, FOV shifted, conjugated, frequency offset corrected.
"""
function water_kdata(c; channel_mixing_matrix, frequency_drift=nothing, gradient_delay_us,
                     do_freq_cor=true, fov_shift_rotation=:siemens)
    kdata = MRSI.read_data(c, channel_mixing_matrix, frequency_drift, ComplexF32)
    MRSI.fov_shift!(kdata, c, gradient_delay_us; fov_shift_rotation)
    kdata = conj.(kdata)
    do_freq_cor && MRSI.frequency_offset_correction!(kdata, c; fill_factor=2, unalias_lipids=true)
    return kdata
end

"""
    circle_onto_grid(kdata, c, grid; gradient_delay_us, do_hamming_filter=true)

The in-plane image of one water circle on the metabolite matrix, (n, n, points, channels):
density compensated with the metabolite's effective area and in-plane window at the
circle's radius (ICE convention, points averaged), transformed at the circle's own
trajectory.

The ICE density compensation weighs each circle by its nominal area, so a scan whose
circles were acquired at radii a little off the nominal ones carries the square of that
ratio in its gain. Scaling by (water radius / metabolite radius)^2 gives the water the
metabolite's gain.
"""
function circle_onto_grid(kdata, c, grid::MetaboliteGrid; gradient_delay_us, do_hamming_filter=true)
    xy = MRSI.read_first_kspace_coordinate_normalized(c) * grid.k_scale
    dcf = grid.area[c[:circle]] * (abs(xy) / grid.radii[c[:circle]])^2 / c[:n_points_on_circle]
    do_hamming_filter && (dcf *= MRSI.hamming_filter(abs(xy)))
    coordinates = MRSI.construct_circle_coordinates(xy, c[:n_points_on_circle],
                                                    MRSI.gradient_delay_to_rad(gradient_delay_us, c))
    dft = ComplexF32.(MRSI.calculate_dft_matrix(coordinates, grid.n_grid))
    image = dft * reshape(kdata .* Float32(dcf), size(kdata, 1), :)
    return reshape(image, grid.n_grid, grid.n_grid, size(kdata)[2:end]...)
end

"""Partition kz of a water circle, 1/mm, the partition FFT's convention."""
circle_kz(c) = (c[:part] - (c[:n_part] ÷ 2 + 1)) / c[:fov_slice]

"""
    water_onto_metabolite_grid(water_scan, met_scan; n_points, channel_mixing_matrix, do_hamming_filter=true,
                               gradient_delay_us=nothing, fov_shift_rotation=:siemens)

The uncombined water ONLINE scan of `water_scan` (MeasurementInfo) on the grid of `met_scan`,
(nx, ny, nz, n_points, channels), reconstructed in ICE mode with the metabolite scan's point
spread function. The frequency offset correction, the trajectory delays and the drift
correction follow the water protocol as `reconstruct_like_ICE` reads them; `gradient_delay_us`
overrides the delays as there.
"""
function water_onto_metabolite_grid(water_scan, met_scan; n_points, channel_mixing_matrix,
                                    do_hamming_filter=true, gradient_delay_us=nothing,
                                    fov_shift_rotation=:siemens)
    water, met = water_scan[:ONLINE], met_scan[:ONLINE]
    grid = MetaboliteGrid(met, water)
    n_points <= water[:n_fid] || error("The water reference has $(water[:n_fid]) FID points, not $n_points")
    n_circles(water) >= maximum(grid.circles) ||
        error("The water reference covers $(n_circles(water)) circles, fewer than the metabolite scan's " *
              "$(maximum(grid.circles)): its point spread function cannot be matched.")
    (water[:n_part] ÷ 2) / water[:fov_slice] >= (met[:n_part] ÷ 2) / met[:fov_slice] ||
        error("The water reference's kz extent is smaller than the metabolite scan's: its point spread " *
              "function cannot be matched.")

    # reconstruct_like_ICE's protocol gating
    if isnothing(gradient_delay_us)
        gradient_delay_us = water_scan[:protocol_do_grad_delay_correction] ? MRSI.GRADIENT_DELAY_US : MRSI.NO_GRADIENT_DELAY_US
    end
    do_freq_cor = water_scan[:protocol_do_freq_offset_correction]
    drift_mode = water_scan[:protocol_do_freq_drift_correction] ? :per_TI : nothing
    prescan_kw = MRSI.handle_prescans(water_scan; frequency_drift_correction=drift_mode, channel_mixing_matrix)
    frequency_drift = get(prescan_kw, :frequency_drift, nothing)

    circles = filter(c -> in_metabolite_extent(grid, c[:circle], circle_kz(c)), MRSI.sort_into_circles(water, 1))
    kz_list = sort(unique(circle_kz.(circles)))
    slot = Dict(kz => i for (i, kz) in enumerate(kz_list))
    n_channels = size(channel_mixing_matrix, 2)
    println("Julia: water reference onto the metabolite grid, $(length(circles)) of its circles in " *
            "$(length(kz_list)) partitions, $n_points FID points, $n_channels channels")

    per_partition = zeros(ComplexF32, grid.n_grid, grid.n_grid, length(kz_list), n_points, n_channels)
    MRSI.accumulate_into_grid!(per_partition, [slot[circle_kz(c)] for c in circles], circles;
                               budget=1024^3, desc="water onto the metabolite grid: ") do c
        kdata = water_kdata(c; channel_mixing_matrix, frequency_drift, gradient_delay_us, do_freq_cor, fov_shift_rotation)
        circle_onto_grid(view(kdata, :, 1:n_points, :), c, grid; gradient_delay_us, do_hamming_filter)
    end

    E = ComplexF32.(z_operator(grid, kz_list, water[:fov_slice]))
    return apply_z_operator(E, per_partition)
end

n_circles(info) = info[:n_frequency] ÷ 2

"""`E` along dimension 3 of `data` (x, y, kz, points, channels), one point and channel at a time."""
function apply_z_operator(E, data)
    nx, ny, _, n_points, n_channels = size(data)
    result = zeros(eltype(data), nx, ny, size(E, 1), n_points, n_channels)
    Et = transpose(E)
    for ch in 1:n_channels, t in 1:n_points
        slab = reshape(view(data, :, :, :, t, ch), nx * ny, :)
        result[:, :, :, t, ch] .= reshape(slab * Et, nx, ny, :)
    end
    return result
end

"""
    resample_fid(fid, dwell_in, dwell_out)

A FID sampled every `dwell_in`, at the times 0, dwell_out, ... up to its last sample, along
the last dimension. Each voxel is demodulated by its own frequency, the slowly varying
envelope interpolated with a not-a-knot cubic spline and remodulated: exact at the samples,
and accurate between them while the envelope changes little over one input dwell.
"""
function resample_fid(fid::AbstractArray{<:Complex}, dwell_in, dwell_out)
    n_in = size(fid, ndims(fid))
    n_out = floor(Int, (n_in - 1) * dwell_in / dwell_out + 1e-9) + 1
    t_in = (0:n_in-1) .* dwell_in
    t_out = (0:n_out-1) .* dwell_out
    S = spline_operator(n_in, n_out, dwell_out / dwell_in)
    flat = reshape(fid, :, n_in)
    out = zeros(eltype(fid), size(flat, 1), n_out)
    k = min(10, n_in - 1)
    for v in axes(flat, 1)
        y = ComplexF64.(view(flat, v, :))
        any(!iszero, y) || continue
        f = angle(sum(view(y, 2:k+1) .* conj.(view(y, 1:k)))) / (2pi * dwell_in)
        envelope = y .* cispi.(-2f .* t_in)
        out[v, :] .= (S * envelope) .* cispi.(2f .* t_out)
    end
    return reshape(out, size(fid)[1:end-1]..., n_out)
end

"""
The (n_out, n_in) matrix of a not-a-knot cubic spline through n_in unit-spaced samples,
evaluated at 0, step, 2 step, ... (in units of the input spacing).
"""
function spline_operator(n_in, n_out, step)
    n_in >= 4 || error("a spline needs four points")
    # second derivatives M = D y: the not-a-knot conditions (third derivative continuous at
    # the second and the second-to-last knot) close the usual tridiagonal system
    A = zeros(n_in, n_in)
    B = zeros(n_in, n_in)
    A[1, 1:3] .= (1, -2, 1)
    A[n_in, n_in-2:n_in] .= (1, -2, 1)
    for i in 2:n_in-1
        A[i, i-1:i+1] .= (1, 4, 1)
        B[i, i-1:i+1] .= (6, -12, 6)
    end
    D = A \ B
    S = zeros(n_out, n_in)
    for (row, t) in enumerate((0:n_out-1) .* step)
        i = clamp(floor(Int, t + 1e-9) + 1, 1, n_in - 1)
        a, d = i - t, t - (i - 1)
        S[row, i] += a
        S[row, i+1] += d
        S[row, :] .+= ((a^3 - a) .* D[i, :] .+ (d^3 - d) .* D[i+1, :]) ./ 6
    end
    return S
end
