# Tests of water_reference_grid.jl: a synthetic object of known density through the metabolite
# reconstruction (MRSI.jl's own steps) and through the water route onto the metabolite grid must
# give the same image, with the water slab thicker, finer and off centre.
#
#     julia --project=<environment with MRSI> tests/test_water_reference_grid.jl
#
# The two geometries are built from MRSI.jl's test.dat (one partition), its circles copied into
# the partitions of a slab; MRSI_TEST_DAT names another copy of it. Skipped when it is missing.
using Test
using MRSI
using LinearAlgebra
using Statistics
include(joinpath(@__DIR__, "..", "water_reference_grid.jl"))

const TEST_DAT = get(ENV, "MRSI_TEST_DAT", joinpath(dirname(pathof(MRSI)), "..", "test", "data", "test.dat"))
const NO_DELAY = [0.0im, 0.0im, 0.0im]

"""A CRT slab from the 2D test.dat: n_grid matrix, n_part partitions over fov_slice mm, its centre
`offset` mm along the normal from the original's."""
function slab_scan(n_grid, n_part, fov_slice; offset=0.0)
    source = MRSI.read_scan_info(TEST_DAT)[:ONLINE]
    info = MRSI.ScanInfo(deepcopy(source.general), deepcopy(source.individual))
    info[:n_frequency] = n_grid
    info[:n_part] = n_part
    info[:fov_slice] = fov_slice
    info[:position] = float.(source[:position]) .+ offset .* normalize(float.(source[:slice_normal]))
    counts = MRSI.calculate_n_circles_per_part(n_part, n_grid ÷ 2)
    headers = MRSI.ScanHeaderVD[]
    for (p, n_circles) in enumerate(counts), h in source[:headers]
        h[:circle] <= n_circles || continue
        dims = copy(h.dims)
        dims[Int(MRSI.SEG)] = p - n_part ÷ 2
        push!(headers, MRSI.ScanHeaderVD(h.mask, h.type, dims, h.ice_param, h.data_position))
    end
    info[:headers] = headers
    info[:radii] = MRSI.read_radii(info)
    return info
end

# ellipsoid with a Gaussian edge, density 1, relative to the metabolite slab's centre (z along the normal)
const SEMI = (70.0, 85.0, 14.0)
const CENTRE = (9.0, -12.0, 3.0)
function object_ft(k)
    q = 2pi * sqrt(sum(abs2, k .* SEMI))
    core = q < 1e-6 ? 4pi / 3 : 4pi * (sin(q) - q * cos(q)) / q^3
    return prod(SEMI) * core * exp(-2pi^2 * 9 * sum(abs2, k)) * cispi(2 * sum(k .* CENTRE))
end

"""The object's k-space signal on circle c, after conj, in a slab whose centre lies `shift` mm along
the normal from the metabolite slab's."""
function circle_signal(c, shift)
    n = c[:n_frequency]
    xy = MRSI.construct_circle_coordinates(MRSI.read_first_kspace_coordinate_normalized(c), c[:n_points_on_circle], 0.0im)
    kz = circle_kz(c)
    data = [object_ft((real(k) * n / c[:fov_readout], imag(k) * n / c[:fov_phase], kz)) * cispi(-2 * kz * shift) for k in xy]
    return reshape(ComplexF32.(data), :, 1, 1)
end

function metabolite_image(met; kspace_window=:hamming)
    image = zeros(ComplexF32, met[:n_frequency], met[:n_frequency], met[:n_part], 1, 1)
    for c in MRSI.sort_into_circles(met, 1)
        kdata = circle_signal(c, 0.0)
        MRSI.density_compensation!(kdata, c; do_hamming_filter=true, do_hamming_filter_z=:ice, ice=true, kspace_window)
        image[:, :, c[:part], :, :] .+= MRSI.fourier_transform(kdata, c, NO_DELAY)
    end
    MRSI.fft_slice_dim!(image)
    return real.(image[:, :, :, 1, 1])
end

function water_image(met, water; kspace_window=:hamming)
    grid = MetaboliteGrid(met, water; kspace_window)
    circles = filter(c -> in_metabolite_extent(grid, c[:circle], circle_kz(c)), MRSI.sort_into_circles(water, 1))
    kz_list = sort(unique(circle_kz.(circles)))
    per_partition = zeros(ComplexF32, grid.n_grid, grid.n_grid, length(kz_list), 1, 1)
    shift = slab_offset_mm(met, water)
    for c in circles
        j = findfirst(==(circle_kz(c)), kz_list)
        per_partition[:, :, j, :, :] .+= circle_onto_grid(circle_signal(c, shift), c, grid; gradient_delay_us=NO_DELAY)
    end
    return real.(apply_z_operator(ComplexF32.(z_operator(grid, kz_list, water[:fov_slice])), per_partition)[:, :, :, 1, 1])
end

"""Centroid along `axis`, in voxels, of the part of `image` above half its maximum."""
function centroid(image, axis)
    weight = image .* (image .> maximum(image) / 2)
    return sum(weight .* reshape(axes(image, axis), ntuple(d -> d == axis ? size(image, axis) : 1, 3))) / sum(weight)
end

if !isfile(TEST_DAT)
    @warn "No $TEST_DAT, the water reference grid tests are skipped"
else
    @testset "water reference with another slab" begin
        met = slab_scan(14, 7, 40.0)
        reference = metabolite_image(met)
        interior = reference .> 0.97 * maximum(reference)
        @test count(interior) > 20
        for offset in (0.0, 5.0, -3.5)
            water = slab_scan(16, 11, 55.0; offset)
            @test slab_offset_mm(met, water) ≈ offset atol = 1e-9
            image = water_image(met, water)
            ratio = image[interior] ./ reference[interior]
            @test median(ratio) ≈ 1 atol = 0.01
            @test maximum(abs.(ratio .- 1)) < 0.03
            # positions in voxels of the metabolite grid
            for axis in 1:3
                @test centroid(image, axis) ≈ centroid(reference, axis) atol = 0.02
            end
        end
        # the protocol states no window here, so the grid takes MRSI.jl's default
        @test MetaboliteGrid(met, slab_scan(16, 11, 55.0)).kspace_window === :hamming
        # the water follows the metabolite reconstruction's k-space window
        for kspace_window in (:hann, :none)
            windowed = metabolite_image(met; kspace_window)
            inside = windowed .> 0.97 * maximum(windowed)
            image = water_image(met, slab_scan(16, 11, 55.0; offset=5.0); kspace_window)
            ratio = image[inside] ./ windowed[inside]
            @test median(ratio) ≈ 1 atol = 0.01
            @test maximum(abs.(ratio .- 1)) < 0.03
        end
        @test_throws ArgumentError MetaboliteGrid(met, slab_scan(16, 11, 55.0); kspace_window=:tukey)
        # a water slab that does not cover the metabolite's kz extent cannot match it
        thin = slab_scan(16, 7, 55.0)
        @test_throws ErrorException water_onto_metabolite_grid(MRSI.MeasurementInfo(Dict(:ONLINE => thin)),
                                                                MRSI.MeasurementInfo(Dict(:ONLINE => met));
                                                                n_points=1, channel_mixing_matrix=I(1))
        # another in-plane centre is not this route
        shifted = slab_scan(16, 11, 55.0)
        shifted[:position] = float.(shifted[:position]) .+ [3.0, 0.0, 0.0]
        @test_throws ErrorException slab_offset_mm(met, shifted)
    end
end

@testset "FID resampling" begin
    dwell_in, dwell_out = 1.12e-3, 0.36e-3
    fid(t) = 3 * cispi(2 * 37.0 * t + 0.4) * exp(-pi * 15.0 * t - (11.0 * t)^2)
    t_in = (0:143) .* dwell_in
    resampled = resample_fid(reshape(fid.(t_in), 1, 1, 1, :), dwell_in, dwell_out)
    @test size(resampled) == (1, 1, 1, 445)
    t_out = (0:444) .* dwell_out
    @test maximum(abs.(vec(resampled) .- fid.(t_out))) < 1e-5 * 3
    # exact at the common sample times
    @test resampled[1, 1, 1, 1] ≈ fid(0.0)
    @test all(iszero, resample_fid(zeros(ComplexF32, 2, 144), dwell_in, dwell_out))
end
