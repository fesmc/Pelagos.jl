# Unit tests for hydrostatic pressure integration.

using Test
using Pelagos.Pressure: compute_pressure, compute_pressure!
using Pelagos.UNESCO: seawater_density

@testset "Hydrostatic pressure" begin

    # Simple 1-column test: uniform T=10, S=35, 4 layers of 250m each.
    # Expected pressure at layer k centre ≈ ρ·g·depth / 1e5 bar.
    @testset "Single ocean column, uniform T/S" begin
        nlon, nlat, nz = 1, 1, 4
        dz     = 250.0         # m per layer
        z_faces = [0.0, -250.0, -500.0, -750.0, -1000.0]
        T      = fill(10.0, nlon, nlat, nz)
        S      = fill(35.0, nlon, nlat, nz)
        mask   = fill(true,  nlon, nlat)

        p = compute_pressure(T, S, z_faces, mask)

        # Pressure at surface-layer centre (125 m depth): ≈ ρ·g·125/1e5
        rho_ref = seawater_density(10.0, 35.0, 0.0)
        g       = 9.81
        p_expect_k1 = rho_ref * g * 125.0 / 1e5   # bar
        @test isapprox(p[1,1,1], p_expect_k1, rtol=1e-2)

        # Pressure must increase monotonically with depth
        @test p[1,1,1] < p[1,1,2] < p[1,1,3] < p[1,1,4]
    end

    # Land cells must have p = 0
    @testset "Land cells have zero pressure" begin
        T    = fill(10.0, 2, 2, 3)
        S    = fill(35.0, 2, 2, 3)
        # [true false; true true] → mask[1,2]=false (land at i=1, j=2)
        mask = [true false; true true]
        z_faces = [0.0, -100.0, -200.0, -300.0]
        p = compute_pressure(T, S, z_faces, mask)
        @test all(p[1, 2, :] .== 0.0)   # land column at (i=1, j=2)
    end

    # In-place and allocating versions give the same result
    @testset "In-place consistency" begin
        T    = randn(4, 4, 6) .* 5.0 .+ 10.0
        S    = randn(4, 4, 6) .* 2.0 .+ 35.0
        mask = fill(true, 4, 4)
        z_faces = collect(range(0.0, stop=-3000.0, length=7))
        p1 = compute_pressure(T, S, z_faces, mask)
        p2 = zeros(4, 4, 6)
        compute_pressure!(p2, T, S, z_faces, mask)
        @test p1 ≈ p2
    end
end
