# Integration tests for Pelagos.jl Phase 5.
#
# Uses the CLIMBER-X PI steady-state restart in climber-x/restart/pi/ as the
# initial condition and as a loose reference for drift checks.
# Skipped automatically in CI when the climber-x directory is absent.

using Test
using Oceananigans: CPU

const CLIMBERX_DIR    = joinpath(@__DIR__, "..", "..", "climber-x")
const PI_RESTART_DIR  = joinpath(CLIMBERX_DIR, "restart", "pi")
const PI_RESTART_FILE = joinpath(PI_RESTART_DIR, "ocn_restart.nc")
const FORCING_FILE    = joinpath(CLIMBERX_DIR, "input", "fake_clim_const_WFDEI-ERA_preind.nc")

@testset "Integration / Phase 5" begin

    if !isfile(PI_RESTART_FILE)
        @warn "PI restart not found at $PI_RESTART_FILE; skipping integration tests."
        @test_skip "Integration tests skipped (no CLIMBER-X restart)"
        return
    end

    @testset "load_climberx_restart" begin
        r = load_climberx_restart(PI_RESTART_DIR)

        @test size(r.T)     == (72, 36, 23)
        @test size(r.S)     == (72, 36, 23)
        @test size(r.u_B)   == (73, 37, 23)
        @test size(r.v_B)   == (73, 37, 23)
        @test size(r.f_ocn) == (72, 36)
        @test size(r.zw)    == (24,)

        ocean_mask = r.f_ocn .> 0.5
        T_ocean    = r.T[:, :, end][ocean_mask]
        S_ocean    = r.S[:, :, end][ocean_mask]
        @test minimum(T_ocean) > -3.0
        @test maximum(T_ocean) < 35.0
        @test minimum(S_ocean) > 20.0
        @test maximum(S_ocean) < 40.0

        # zw ordered bottom→top (most negative first), surface at 0
        @test issorted(r.zw)
        @test r.zw[end] ≈ 0.0 atol=0.1
    end

    @testset "bgrid → cgrid interpolation" begin
        r   = load_climberx_restart(PI_RESTART_DIR)
        u_C = bgrid_u_to_cgrid(r.u_B)
        v_C = bgrid_v_to_cgrid(r.v_B)

        @test size(u_C) == (73, 36, 23)   # Face in λ, Center in φ
        @test size(v_C) == (72, 37, 23)   # Center in λ, Face in φ

        # Averaging cannot exceed original range
        @test minimum(u_C) >= minimum(r.u_B) - 1e-10
        @test maximum(u_C) <= maximum(r.u_B) + 1e-10
    end

    @testset "load_climberx_forcing" begin
        if !isfile(FORCING_FILE)
            @test_skip "Forcing file not found; skipping"
        else
            f = load_climberx_forcing(FORCING_FILE)
            @test size(f.taux) == (72, 36, 12)
            @test size(f.tauy) == (72, 36, 12)
            @test maximum(abs.(f.taux)) < 5.0   # physically plausible (< 5 N m⁻²)
            @test maximum(abs.(f.tauy)) < 5.0
        end
    end

    @testset "build_climberx_grid" begin
        grid = build_climberx_grid(CPU(); restart_file=PI_RESTART_FILE)
        @test grid.Nx == 72
        @test grid.Ny == 36
        @test grid.Nz == 23
    end

end
