# Tests for the tracer diffusion closures.

using Test
using Pelagos.Diffusion: bryan_lewis_kappa, gm_redi_closure, diapycnal_closure
using Pelagos.Parameters: K_BG, K_DEEP, K_DEEP_ZREF

@testset "Tracer diffusion" begin

    @testset "Bryan-Lewis profile limits" begin
        # κ(z) = K_BG + (K_DEEP - K_BG) * (2/π) * atan(|z| / K_DEEP_ZREF)
        # At z=0:       κ = K_BG  (surface minimum)
        # At z=-z_ref:  κ = (K_BG + K_DEEP) / 2  (half-saturation depth)
        # At z→-∞:      κ → K_DEEP  (deep asymptote)

        @test bryan_lewis_kappa(0.0) ≈ K_BG atol=1e-15

        kappa_half = bryan_lewis_kappa(-K_DEEP_ZREF)
        @test isapprox(kappa_half, (K_BG + K_DEEP) / 2, rtol=1e-10)

        kappa_deep = bryan_lewis_kappa(-1e6)
        @test isapprox(kappa_deep, K_DEEP, rtol=1e-4)

        # Monotonically increasing with depth
        depths = [0.0, -200.0, -500.0, -1000.0, -2000.0, -5000.0]
        kappas = [bryan_lewis_kappa(z) for z in depths]
        @test issorted(kappas)

        # Always in [K_BG, K_DEEP]
        for z in depths
            k = bryan_lewis_kappa(z)
            @test k >= K_BG - 1e-20
            @test k <= K_DEEP + 1e-20
        end
    end

    @testset "Bryan-Lewis: positive diffusivity" begin
        for z in [0.0, -100.0, -1000.0, -4000.0]
            @test bryan_lewis_kappa(z) > 0.0
        end
    end
end
