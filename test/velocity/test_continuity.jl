# Tests for the continuity equation and vertical velocity diagnosis.

using Test
using Pelagos.Continuity: diagnose_w, diagnose_w!

@testset "Vertical velocity from continuity" begin

    @testset "Non-divergent horizontal flow → w = 0 everywhere" begin
        # Uniform u = const, v = 0 → ∂u/∂x = ∂v/∂y = 0 → w = 0
        nlon, nlat, nz = 4, 4, 3
        u    = fill(0.1, nlon, nlat, nz)
        v    = fill(0.0, nlon, nlat, nz)
        cos_lat = fill(1.0, nlat)       # pretend equatorial (cos φ = 1)
        dx   = fill(1e5, nlat)
        dy   = 1e5
        dz   = [100.0, 200.0, 500.0]
        mask = fill(true, nlon, nlat)

        w = diagnose_w(u, v, dx, dy, cos_lat, dz, mask)

        # With periodic longitude and uniform u, ∂u/∂x = 0 everywhere
        @test maximum(abs.(w)) < 1e-10
    end

    @testset "Land cells produce zero w" begin
        nlon, nlat, nz = 3, 3, 2
        u    = zeros(nlon, nlat, nz)
        v    = zeros(nlon, nlat, nz)
        cos_lat = fill(1.0, nlat)
        dx   = fill(5e4, nlat)
        dy   = 5e4
        dz   = [50.0, 100.0]
        mask = fill(true, nlon, nlat)
        mask[2, 2] = false

        w = diagnose_w(u, v, dx, dy, cos_lat, dz, mask; debug=false)
        @test all(w[2, 2, :] .== 0.0)
    end

    @testset "Bottom w = 0" begin
        nlon, nlat, nz = 3, 3, 4
        u    = randn(nlon, nlat, nz) .* 0.01
        v    = randn(nlon, nlat, nz) .* 0.01
        cos_lat = fill(1.0, nlat)
        dx   = fill(5e5, nlat)
        dy   = 5e5
        dz   = [50.0, 100.0, 200.0, 500.0]
        mask = fill(true, nlon, nlat)

        w = diagnose_w(u, v, dx, dy, cos_lat, dz, mask)
        # Bottom face (k = nz+1) must be exactly zero
        @test all(w[:, :, end] .== 0.0)
    end

    @testset "In-place and allocating versions agree" begin
        nlon, nlat, nz = 3, 3, 3
        u    = randn(nlon, nlat, nz) .* 0.01
        v    = randn(nlon, nlat, nz) .* 0.01
        cos_lat = fill(cosd(30.0), nlat)
        dx   = fill(4e5, nlat)
        dy   = 4e5
        dz   = [100.0, 200.0, 400.0]
        mask = fill(true, nlon, nlat)

        w1 = diagnose_w(u, v, dx, dy, cos_lat, dz, mask)
        w2 = zeros(nlon, nlat, nz+1)
        diagnose_w!(w2, u, v, dx, dy, cos_lat, dz, mask)
        @test w1 ≈ w2
    end
end
