# Tests for the tracer diffusion closures.

using Test
using Pelagos.Diffusion: bryan_lewis_kappa, gm_redi_closure, diapycnal_closure
using Pelagos.Parameters: K_BG, K_DEEP

@testset "Tracer diffusion" begin

    @testset "Bryan-Lewis profile limits" begin
        # At z=0, the arctan argument is -α*z_ref (large negative) → κ < K_BG.
        # At z = -z_ref, arctan(0) = 0 → κ = K_BG exactly.
        # At z → -∞, arctan → π/2 → κ → K_DEEP.
        kappa_surface = bryan_lewis_kappa(0.0)
        @test kappa_surface > 0.0
        @test kappa_surface < K_DEEP

        # At z_ref depth, κ should equal K_BG exactly
        using Pelagos.Parameters: K_DEEP_ZREF
        kappa_zref = bryan_lewis_kappa(-K_DEEP_ZREF)
        @test isapprox(kappa_zref, K_BG, atol=1e-12)

        # At great depth (z → -∞), κ → K_DEEP; need z >> z_ref for close convergence
        kappa_deep = bryan_lewis_kappa(-100000.0)
        @test isapprox(kappa_deep, K_DEEP, rtol=0.02)

        # κ must be monotonically increasing with depth
        depths = [0.0, -500.0, -1000.0, -2000.0, -4000.0, -6000.0]
        kappas = [bryan_lewis_kappa(z) for z in depths]
        @test issorted(kappas)
    end

    @testset "Bryan-Lewis: positive diffusivity" begin
        for z in [0.0, -100.0, -1000.0, -4000.0]
            @test bryan_lewis_kappa(z) > 0.0
        end
    end
end
