# Unit tests for hydrostatic pressure integration on an Oceananigans grid.
#
# Phase-5 API: compute_pressure!(p::Field{C,C,C}, T::Field, S::Field, grid).
# The plain-array entry point of the original implementation has been removed
# and tests now construct a small ImmersedBoundaryGrid + Fields directly.

using Test
using Oceananigans
using Oceananigans.Fields: Field, Center, interior, set!
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom

using Pelagos.Pressure: compute_pressure!
using Pelagos.UNESCO: seawater_density

# Build a tiny LatitudeLongitudeGrid wrapped in an ImmersedBoundaryGrid.
# Halos are sized to fit small test grids (production code uses halo=(4,4,4) for
# higher-order stencils, but pressure integration only needs minimal halos).
function _make_test_grid(Nx, Ny, z_faces; bottom_height=nothing)
    ug = LatitudeLongitudeGrid(CPU();
        size      = (Nx, Ny, length(z_faces)-1),
        longitude = (-180, 180),
        latitude  = (-30.0, 30.0),
        z         = z_faces,
        topology  = (Periodic, Bounded, Bounded),
        halo      = (1, 1, 1),
    )
    bh = isnothing(bottom_height) ? fill(z_faces[1], Nx, Ny) : bottom_height
    return ImmersedBoundaryGrid(ug, GridFittedBottom(bh))
end

@testset "Hydrostatic pressure" begin

    # 4 layers, 250 m each.  Surface k=Nz=4 in Oceananigans convention.
    @testset "Single ocean column, uniform T/S" begin
        z_faces = [-1000.0, -750.0, -500.0, -250.0, 0.0]   # bottom→top
        grid    = _make_test_grid(1, 1, z_faces)

        T = Field{Center, Center, Center}(grid)
        S = Field{Center, Center, Center}(grid)
        p = Field{Center, Center, Center}(grid)
        set!(T, 10.0)
        set!(S, 35.0)

        compute_pressure!(p, T, S, grid)
        pi = interior(p)                # (1, 1, 4)

        # Surface-cell centre is at depth 125 m (Nz=4 layer is the top 250 m,
        # center 125 m below z=0).  ρ·g·125/1e5 in bar.
        ρ_ref = seawater_density(10.0, 35.0, 0.0)
        p_expect_top = ρ_ref * 9.81 * 125.0 / 1e5
        @test isapprox(pi[1, 1, 4], p_expect_top, rtol=1e-2)

        # Pressure increases monotonically with depth (k=1 deepest, k=Nz top)
        @test pi[1, 1, 1] > pi[1, 1, 2] > pi[1, 1, 3] > pi[1, 1, 4]
    end

    # In the new API, fully-immersed (land) columns are filled from an adjacent
    # ocean column to suppress spurious horizontal pressure gradients at coasts.
    # Verify the inheritance instead of the old "land = 0" semantics.
    @testset "Land columns inherit pressure from ocean neighbour" begin
        z_faces = [-300.0, -200.0, -100.0, 0.0]
        Nx, Ny = 2, 2
        # Make column (1, 2) fully land (bottom_height = 0 seals the column)
        bh = fill(z_faces[1], Nx, Ny)
        bh[1, 2] = 0.0
        grid = _make_test_grid(Nx, Ny, z_faces; bottom_height=bh)

        T = Field{Center, Center, Center}(grid)
        S = Field{Center, Center, Center}(grid)
        p = Field{Center, Center, Center}(grid)
        set!(T, 10.0)
        set!(S, 35.0)
        compute_pressure!(p, T, S, grid)
        pi = interior(p)

        # Land column must be finite and equal an ocean neighbour's profile,
        # not zero — that's the new "ghost-fill" design.
        for k in 1:size(pi, 3)
            @test isfinite(pi[1, 2, k])
            ocean_neighbours = (pi[2, 2, k], pi[1, 1, k])   # i+1 then j-1
            @test minimum(abs.(pi[1, 2, k] .- ocean_neighbours)) < 1e-10
        end
    end

    # Determinism: repeated calls give identical results.
    @testset "Repeated calls produce the same pressure field" begin
        z_faces = collect(range(-3000.0, 0.0; length=7))
        grid = _make_test_grid(4, 4, z_faces)
        T = Field{Center, Center, Center}(grid)
        S = Field{Center, Center, Center}(grid)
        p1 = Field{Center, Center, Center}(grid)
        p2 = Field{Center, Center, Center}(grid)
        set!(T, (x, y, z) -> 10.0 + 0.001*z)
        set!(S, (x, y, z) -> 35.0 - 0.0005*z)

        compute_pressure!(p1, T, S, grid)
        compute_pressure!(p2, T, S, grid)
        @test interior(p1) ≈ interior(p2)
    end
end
