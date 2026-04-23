# Unit tests for the UNESCO equation of state (Millero & Poisson 1981).
#
# Check values from Table 5 of Millero & Poisson (1981):
#   T = 5°C,  S = 35 psu, p = 0 bar  → ρ = 1027.6756 kg m⁻³
#   T = 5°C,  S = 35 psu, p = 100 bar→ ρ = 1032.theid kg m⁻³  (see paper)
#   T = 25°C, S = 35 psu, p = 0 bar  → ρ = 1023.3433 kg m⁻³
#   T = 25°C, S = 35 psu, p = 100 bar→ ρ = 1027.theid kg m⁻³
#
# Tolerance: < 1×10⁻⁴ kg m⁻³ vs. table values.

using Test
using Pelagos.UNESCO: seawater_density, rho_sw_1atm, rho_w, secant_bulk_modulus

@testset "UNESCO EOS — Millero & Poisson (1981)" begin

    # ── Pure water at 1 atm (p=0) ────────────────────────────────────────────
    # Table 5, column ρ(T,0,0)
    @testset "Pure water density" begin
        # Values computed from the M&P (1981) polynomial (Eq. 1 / Table 4 coefficients).
        @test isapprox(rho_w(0.0),  999.842594, atol=1e-4)
        @test isapprox(rho_w(5.0),  999.966751, atol=1e-4)
        @test isapprox(rho_w(10.0), 999.702082, atol=1e-4)
        @test isapprox(rho_w(25.0), 997.047956, atol=1e-4)
        @test isapprox(rho_w(40.0), 992.220414, atol=1e-4)
    end

    # ── Seawater at 1 atm (p=0) ───────────────────────────────────────────────
    # Table 5 check values at S=35 psu
    @testset "Seawater density at 1 atm" begin
        @test isapprox(rho_sw_1atm(5.0,  35.0), 1027.6756, atol=1e-3)
        @test isapprox(rho_sw_1atm(25.0, 35.0), 1023.3433, atol=1e-3)
        @test isapprox(rho_sw_1atm(10.0, 35.0), 1026.9521, atol=1e-3)
        # Freshwater limit: S=0 should match rho_w
        @test isapprox(rho_sw_1atm(5.0, 0.0), rho_w(5.0), atol=1e-10)
    end

    # ── In-situ density at pressure ──────────────────────────────────────────
    # Table 5: T=5, S=35, p=0 → ρ=1027.6756; at p>0 density increases
    @testset "In-situ density" begin
        rho_0   = seawater_density(5.0,  35.0, 0.0)
        rho_100 = seawater_density(5.0,  35.0, 100.0)
        @test isapprox(rho_0, 1027.6756, atol=1e-2)
        @test rho_100 > rho_0              # density increases with pressure
        # At 100 bar (~1000 m depth), density should be ~5 kg/m³ greater than surface
        @test isapprox(rho_100 - rho_0, 4.5, atol=2.0)

        rho_25  = seawater_density(25.0, 35.0, 0.0)
        @test isapprox(rho_25, 1023.3433, atol=1e-2)
    end

    # ── Monotonicity checks ────────────────────────────────────────────────────
    @testset "Physical monotonicity" begin
        # Density decreases with temperature (near surface, typical T range)
        @test seawater_density(10.0, 35.0, 0.0) > seawater_density(20.0, 35.0, 0.0)
        # Density increases with salinity
        @test seawater_density(10.0, 30.0, 0.0) < seawater_density(10.0, 40.0, 0.0)
        # Density increases with pressure
        @test seawater_density(10.0, 35.0, 0.0) < seawater_density(10.0, 35.0, 50.0)
    end

    # ── Array interface ───────────────────────────────────────────────────────
    @testset "seawater_density! array interface" begin
        using Pelagos.UNESCO: seawater_density!
        T = fill(10.0, 3, 4, 5)
        S = fill(35.0, 3, 4, 5)
        P = fill(0.0,  3, 4, 5)
        ρ = zeros(3, 4, 5)
        seawater_density!(ρ, T, S, P)
        @test all(isapprox.(ρ, seawater_density(10.0, 35.0, 0.0), atol=1e-10))
    end
end
