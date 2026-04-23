# Tests for the frictional-geostrophic baroclinic velocity solver.

using Test
using Pelagos.Baroclinic: solve_baroclinic!, coriolis_parameter
using Pelagos.Parameters: R_BC, F_MIN, RHO_0, OMEGA

@testset "Baroclinic velocity solver" begin

    @testset "Coriolis parameter with floor" begin
        # At 45°N: f = 2Ω·sin(45°)
        f45 = coriolis_parameter(45.0)
        @test isapprox(f45, 2*OMEGA*sind(45.0), rtol=1e-10)

        # Equatorial floor: |f| must be ≥ F_MIN
        f_eq = coriolis_parameter(0.0)
        @test abs(f_eq) ≥ F_MIN

        f_near_eq = coriolis_parameter(0.01)
        @test abs(f_near_eq) ≥ F_MIN
    end

    @testset "Analytical solution: no pressure gradient, only wind" begin
        # With no pressure gradient, the balance is:
        #   -f·v + r·u = τˣ/(ρ·h)
        #    f·u + r·v = τʸ/(ρ·h)
        # At 45°N with τˣ = 0.1 N m⁻², τʸ = 0, h = 50 m:
        nlon, nlat, nz = 3, 3, 1
        u      = zeros(nlon, nlat, nz)
        v      = zeros(nlon, nlat, nz)
        p      = zeros(nlon, nlat, nz)   # no pressure gradient
        tau_x  = fill(0.1, nlon, nlat)   # N m⁻²
        tau_y  = fill(0.0, nlon, nlat)
        f      = [coriolis_parameter(45.0) for _ in 1:nlat]
        dx     = fill(5.0e5, nlat)   # 500 km
        dy     = 5.0e5
        dz     = [50.0]
        mask   = fill(true, nlon, nlat)

        solve_baroclinic!(u, v, p, tau_x, tau_y, f, dx, dy, dz, mask)

        fj  = f[2]
        h   = dz[1]
        Fx  = tau_x[2,2] / (RHO_0 * h)
        Fy  = 0.0
        r   = R_BC
        u_exp = (r*Fx + fj*Fy) / (r^2 + fj^2)
        v_exp = (r*Fy - fj*Fx) / (r^2 + fj^2)

        @test isapprox(u[2,2,1], u_exp, rtol=1e-8)
        @test isapprox(v[2,2,1], v_exp, rtol=1e-8)
    end

    @testset "Land cells produce zero velocity" begin
        nlon, nlat, nz = 4, 4, 2
        u      = fill(999.0, nlon, nlat, nz)
        v      = fill(999.0, nlon, nlat, nz)
        p      = zeros(nlon, nlat, nz)
        tau_x  = zeros(nlon, nlat)
        tau_y  = zeros(nlon, nlat)
        f      = [coriolis_parameter(30.0) for _ in 1:nlat]
        dx     = fill(5e5, nlat)
        dy     = 5e5
        dz     = [50.0, 100.0]
        mask   = fill(true, nlon, nlat)
        mask[2, 2] = false   # one land cell

        solve_baroclinic!(u, v, p, tau_x, tau_y, f, dx, dy, dz, mask)
        @test u[2, 2, 1] == 0.0
        @test v[2, 2, 1] == 0.0
        @test u[2, 2, 2] == 0.0
    end

    @testset "Below surface, wind forcing has no effect" begin
        nlon, nlat, nz = 2, 2, 3
        u      = zeros(nlon, nlat, nz)
        v      = zeros(nlon, nlat, nz)
        p      = zeros(nlon, nlat, nz)
        tau_x  = fill(0.2, nlon, nlat)
        tau_y  = fill(0.1, nlon, nlat)
        f      = [coriolis_parameter(60.0) for _ in 1:nlat]
        dx     = fill(3e5, nlat)
        dy     = 3e5
        dz     = [50.0, 100.0, 200.0]
        mask   = fill(true, nlon, nlat)

        solve_baroclinic!(u, v, p, tau_x, tau_y, f, dx, dy, dz, mask)

        # At k=2 and k=3 the only forcing is the pressure gradient (zero here)
        @test u[1,1,2] == 0.0
        @test u[1,1,3] == 0.0
    end
end
