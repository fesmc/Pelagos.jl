# Tests for the barotropic streamfunction solver.

using Test
using SparseArrays
using Pelagos.Barotropic: build_barotropic_solver, solve_barotropic!, compute_rbt
using Pelagos.Baroclinic: coriolis_parameter

@testset "Barotropic streamfunction solver" begin

    @testset "compute_rbt — enhancement near boundaries" begin
        # All-ocean grid so only the shallow-water criterion fires.
        nlon, nlat = 5, 5
        ocean_mask = fill(true, nlon, nlat)

        H    = fill(3000.0, nlon, nlat)
        H[1, :] .= 300.0   # shallow western column (H < 500 m threshold)
        dx   = fill(5e5, nlat)
        dy   = 5e5

        r_bt = compute_rbt(ocean_mask, H, dx, dy)

        # Column i=1 is shallow → enhanced friction; column i=3 (deep interior) is not.
        for j in 1:nlat
            @test r_bt[1, j] > r_bt[3, j]
        end
    end

    @testset "Solver produces finite solution on simple domain" begin
        nlon, nlat = 8, 6
        ocean_mask = fill(true, nlon, nlat)
        # Add a continent along east boundary
        ocean_mask[7:8, :] .= false

        H    = fill(3000.0, nlon, nlat)
        H[ocean_mask .== false] .= 0.0
        f    = [coriolis_parameter(Float64(j * 10 - 30)) for j in 1:nlat]
        dx   = fill(5e5, nlat)
        dy   = 5e5
        r_bt = compute_rbt(ocean_mask, H, dx, dy)

        solver = build_barotropic_solver(ocean_mask, f, H, dx, dy, r_bt)

        tau_x  = fill(0.1, nlon, nlat)
        tau_y  = fill(0.0, nlon, nlat)

        psi = solve_barotropic!(solver, tau_x, tau_y, H, f, dx, dy, ocean_mask)

        @test size(psi) == (nlon, nlat)
        @test all(isfinite.(psi))
        # Land cells must have ψ = 0
        @test all(psi[.!ocean_mask] .== 0.0)
    end

    @testset "Zero wind stress → near-zero streamfunction" begin
        nlon, nlat = 6, 4
        ocean_mask = fill(true, nlon, nlat)
        H    = fill(4000.0, nlon, nlat)
        f    = [coriolis_parameter(Float64(j * 15)) for j in 1:nlat]
        dx   = fill(4e5, nlat)
        dy   = 4e5
        r_bt = compute_rbt(ocean_mask, H, dx, dy)

        solver = build_barotropic_solver(ocean_mask, f, H, dx, dy, r_bt)

        tau_x = fill(0.0, nlon, nlat)
        tau_y = fill(0.0, nlon, nlat)
        psi   = solve_barotropic!(solver, tau_x, tau_y, H, f, dx, dy, ocean_mask)

        @test maximum(abs.(psi)) < 1e-6
    end
end
