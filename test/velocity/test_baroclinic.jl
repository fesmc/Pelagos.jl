# Tests for the frictional-geostrophic baroclinic velocity solver.
#
# Phase-5 API: solve_baroclinic!(u::Field{Face,C,C}, v::Field{C,Face,C},
#                                p::Field{C,C,C}, tau_x::Matrix, tau_y::Matrix,
#                                grid::ImmersedBoundaryGrid).

using Test
using Oceananigans
using Oceananigans.Fields: Field, Center, Face, interior, set!
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom

using Pelagos.Baroclinic: solve_baroclinic!, coriolis_parameter
using Pelagos.Parameters: R_BC, F_MIN, RHO_0, OMEGA

# Build a small test ImmersedBoundaryGrid centred at a given mid-latitude.
# `lat_centre` controls the latitude of the band so we can target a known f.
function _make_baroclinic_grid(Nx, Ny, z_faces; lat_centre=45.0, bottom_height=nothing)
    half = 2.5     # ±2.5° band — keeps Δlat small so f is ~uniform across the grid
    ug = LatitudeLongitudeGrid(CPU();
        size      = (Nx, Ny, length(z_faces)-1),
        longitude = (-180, 180),
        latitude  = (lat_centre - half, lat_centre + half),
        z         = z_faces,
        topology  = (Periodic, Bounded, Bounded),
        halo      = (1, 1, 1),
    )
    bh = isnothing(bottom_height) ? fill(z_faces[1], Nx, Ny) : bottom_height
    return ImmersedBoundaryGrid(ug, GridFittedBottom(bh))
end

@testset "Baroclinic velocity solver" begin

    @testset "Coriolis parameter with floor" begin
        f45 = coriolis_parameter(45.0)
        @test isapprox(f45, 2*OMEGA*sind(45.0), rtol=1e-10)

        @test abs(coriolis_parameter(0.0))  ≥ F_MIN
        @test abs(coriolis_parameter(0.01)) ≥ F_MIN
    end

    # No pressure gradient → wind alone drives the surface flow according to
    #   u = (r·Fx + f·Fy)/(r² + f²)
    #   v = (r·Fy − f·Fx)/(r² + f²)
    # with Fx = τx/(ρ·h), Fy = 0, h = top-layer thickness.
    @testset "Analytical solution: no pressure gradient, only wind" begin
        z_faces = [-200.0, -150.0, -100.0, -50.0, 0.0]   # 4 layers, 50 m each
        Nx, Ny = 4, 4
        grid = _make_baroclinic_grid(Nx, Ny, z_faces; lat_centre=45.0)

        u = Field{Face,   Center, Center}(grid)
        v = Field{Center, Face,   Center}(grid)
        p = Field{Center, Center, Center}(grid)
        # p ≡ 0 → ∂p/∂x = ∂p/∂y = 0
        tau_x = fill(0.1, Nx, Ny)
        tau_y = fill(0.0, Nx, Ny)

        solve_baroclinic!(u, v, p, tau_x, tau_y, grid)
        ui = interior(u); vi = interior(v)

        # Solver evaluates Coriolis at the U-face latitude (T-row centre) and
        # at the V-face latitude (T-row interface).  Compute expected values
        # using the same grid-point latitudes for an exact match.
        i, j = 2, 2
        Nz   = grid.Nz
        ug   = grid.underlying_grid
        f_u  = coriolis_parameter(ug.φᵃᶜᵃ[j])
        f_v  = coriolis_parameter(ug.φᵃᶠᵃ[j])
        r    = R_BC
        h    = ug.z.Δᵃᵃᶜ[Nz]
        Fx   = 0.1 / (RHO_0 * h)
        u_exp = (r*Fx)        / (r^2 + f_u^2)
        v_exp = (-f_v*Fx)     / (r^2 + f_v^2)

        @test isapprox(ui[i, j, Nz], u_exp, rtol=1e-6)
        @test isapprox(vi[i, j, Nz], v_exp, rtol=1e-6)
    end

    # An immersed (land) cell has its surrounding faces marked peripheral, so
    # the solver must leave u, v at those faces equal to zero.
    @testset "Land cells produce zero velocity" begin
        z_faces = [-100.0, -50.0, 0.0]    # 2 layers
        Nx, Ny = 4, 4
        bh = fill(z_faces[1], Nx, Ny)
        bh[2, 2] = 0.0                    # column (2, 2) is land (sealed top to bottom)
        grid = _make_baroclinic_grid(Nx, Ny, z_faces; lat_centre=30.0, bottom_height=bh)

        u = Field{Face,   Center, Center}(grid)
        v = Field{Center, Face,   Center}(grid)
        p = Field{Center, Center, Center}(grid)
        # Pre-fill with a sentinel — solver must overwrite or zero faces of land cell.
        fill!(interior(u), 999.0); fill!(interior(v), 999.0)

        solve_baroclinic!(u, v, p, zeros(Nx, Ny), zeros(Nx, Ny), grid)
        ui = interior(u); vi = interior(v)

        # Faces touching the land cell (i=2, j=2) must be zero at every level.
        for k in 1:grid.Nz
            @test ui[2, 2, k] == 0.0    # west face of land cell
            @test ui[3, 2, k] == 0.0    # east face of land cell
            @test vi[2, 2, k] == 0.0    # south face of land cell
            @test vi[2, 3, k] == 0.0    # north face of land cell
        end
    end

    # Wind stress is added only at k = Nz (surface layer).  With no pressure
    # gradient anywhere, deeper layers must have u = v = 0 exactly.
    @testset "Below surface, wind forcing has no effect" begin
        z_faces = [-350.0, -150.0, -50.0, 0.0]   # 3 layers
        Nx, Ny = 4, 4
        grid = _make_baroclinic_grid(Nx, Ny, z_faces; lat_centre=60.0)

        u = Field{Face,   Center, Center}(grid)
        v = Field{Center, Face,   Center}(grid)
        p = Field{Center, Center, Center}(grid)
        tau_x = fill(0.2, Nx, Ny)
        tau_y = fill(0.1, Nx, Ny)

        solve_baroclinic!(u, v, p, tau_x, tau_y, grid)
        ui = interior(u); vi = interior(v)
        Nz = grid.Nz

        for k in 1:(Nz-1), j in 1:Ny, i in 1:Nx
            @test ui[i, j, k] == 0.0
        end
        for k in 1:(Nz-1), j in 2:Ny, i in 1:Nx
            @test vi[i, j, k] == 0.0
        end
    end
end
