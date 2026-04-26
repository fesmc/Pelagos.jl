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

    # Order-of-magnitude Sverdrup balance check.
    #
    # Interior steady-state vorticity balance (uniform H, no friction, β-plane):
    #     β ∂ψ/∂x = − curl(τ) / ρ₀
    # Integrating from a closed eastern boundary (ψ=0):
    #     ψ ∼ (curl(τ) / (β · ρ₀)) · L_x
    #
    # For τ₀ ≈ 0.1 N m⁻², β ≈ 2×10⁻¹¹ (m s)⁻¹, basin width L_x ≈ 4×10⁶ m, the
    # Sverdrup transport is on the order of single-digit Sv.  A missing 1/ρ₀
    # factor on the RHS makes ψ a factor ρ₀ ≈ 1000 too large — easily caught
    # by a 100-Sv ceiling.
    @testset "Sverdrup-balance scaling (order-of-magnitude)" begin
        nlon, nlat = 8, 6
        ocean_mask = fill(true, nlon, nlat)
        ocean_mask[end, :] .= false        # closed eastern wall
        H = fill(4000.0, nlon, nlat)
        H[.!ocean_mask] .= 0.0

        # β-plane: f varies linearly in y, centred near 30°N (subtropical gyre).
        f0  = coriolis_parameter(30.0)
        β   = 2 * 7.2921e-5 * cosd(30.0) / 6.371e6
        dy  = 5.5e5                          # m, ~5° at mid-latitude
        f   = [f0 + β * (j - (nlat+1)/2) * dy for j in 1:nlat]
        dx  = fill(4.8e5, nlat)              # m, ~5° × cos(30°) at mid-latitude

        r_bt   = compute_rbt(ocean_mask, H, dx, dy)
        solver = build_barotropic_solver(ocean_mask, f, H, dx, dy, r_bt)

        # Anticyclonic-gyre forcing: τx peaks in the middle of the band, falls
        # to zero at the north and south walls, so curl(τ) ≠ 0 in the interior.
        τ0    = 0.1                          # N m⁻²
        tau_x = [τ0 * sin(π * (j-1) / (nlat-1)) for _ in 1:nlon, j in 1:nlat]
        tau_y = zeros(nlon, nlat)

        ψ = solve_barotropic!(solver, tau_x, tau_y, H, f, dx, dy, ocean_mask)

        @test all(isfinite.(ψ))
        ψmax = maximum(abs.(ψ))

        # Hard ceiling: a missing /ρ₀ would put ψmax in the 10⁹–10¹⁰ m³/s range
        # (1000–10 000 Sv).  A correctly scaled solver at 5° × 5° resolution
        # gives |ψ| ≲ a few Sv (friction-dominated Stommel-like response, well
        # below the inviscid Sverdrup estimate).
        @test ψmax < 1e8                     # i.e. < 100 Sv
        # Sanity floor: solver must respond to the forcing, not just return ~0.
        @test ψmax > 1e5                     # i.e. > 0.1 Sv

        # Diagnostic only: compare with the inviscid Sverdrup transport.  At
        # this resolution friction is non-negligible, so ψmax sits well below
        # this scale; we report it but don't assert a tight bound.
        curl_τ = π * τ0 / ((nlat-1) * dy)
        Lx     = (nlon-1) * dx[1]
        ψ_sv   = curl_τ / (β * 1000.0) * Lx
        @test ψmax < 5 * ψ_sv                # never larger than ~Sverdrup
    end
end
