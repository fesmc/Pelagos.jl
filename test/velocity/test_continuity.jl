# Tests for the continuity equation and vertical velocity diagnosis.
#
# Phase-5 API: diagnose_w!(w::Field{C,C,Face}, u::Field{Face,C,C},
#                          v::Field{C,Face,C}, grid::ImmersedBoundaryGrid).

using Test
using Oceananigans
using Oceananigans.Fields: Field, Center, Face, interior, set!
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom

using Pelagos.Continuity: diagnose_w!

function _make_continuity_grid(Nx, Ny, z_faces; bottom_height=nothing)
    ug = LatitudeLongitudeGrid(CPU();
        size      = (Nx, Ny, length(z_faces)-1),
        longitude = (-180, 180),
        latitude  = (-2.5, 2.5),
        z         = z_faces,
        topology  = (Periodic, Bounded, Bounded),
        halo      = (1, 1, 1),
    )
    bh = isnothing(bottom_height) ? fill(z_faces[1], Nx, Ny) : bottom_height
    return ImmersedBoundaryGrid(ug, GridFittedBottom(bh))
end

@testset "Vertical velocity from continuity" begin

    # Uniform u, v=0 on a Periodic-x grid → ∂u/∂x = 0 → w stays at zero.
    @testset "Non-divergent horizontal flow → w = 0 everywhere" begin
        z_faces = [-800.0, -300.0, -100.0, 0.0]
        grid = _make_continuity_grid(4, 4, z_faces)
        u = Field{Face,   Center, Center}(grid)
        v = Field{Center, Face,   Center}(grid)
        w = Field{Center, Center, Face  }(grid)
        set!(u, (x, y, z) -> 0.1)        # uniform 0.1 m/s
        fill_halo_regions!(u)
        fill_halo_regions!(v)

        diagnose_w!(w, u, v, grid)
        @test maximum(abs.(interior(w))) < 1e-10
    end

    # Land columns must have w identically zero — diagnose_w! skips inactive
    # nodes via inactive_node().
    @testset "Land cells produce zero w" begin
        z_faces = [-150.0, -50.0, 0.0]
        Nx, Ny = 4, 4
        bh = fill(z_faces[1], Nx, Ny)
        bh[2, 2] = 0.0                   # column (2,2) sealed → land
        grid = _make_continuity_grid(Nx, Ny, z_faces; bottom_height=bh)

        u = Field{Face,   Center, Center}(grid)
        v = Field{Center, Face,   Center}(grid)
        w = Field{Center, Center, Face  }(grid)
        set!(u, (x, y, z) -> 0.05 * sin(x))
        set!(v, (x, y, z) -> 0.03 * cos(y))
        fill_halo_regions!(u); fill_halo_regions!(v)

        diagnose_w!(w, u, v, grid)
        wi = interior(w)
        @test all(wi[2, 2, :] .== 0.0)
    end

    # Bottom face (k=1) is always zero — that's the rigid-bottom integration
    # constant in diagnose_w!.
    @testset "Bottom w = 0" begin
        z_faces = [-800.0, -400.0, -200.0, -50.0, 0.0]
        grid = _make_continuity_grid(4, 4, z_faces)
        u = Field{Face,   Center, Center}(grid)
        v = Field{Center, Face,   Center}(grid)
        w = Field{Center, Center, Face  }(grid)
        set!(u, (x, y, z) -> 0.01 * randn())
        set!(v, (x, y, z) -> 0.01 * randn())
        fill_halo_regions!(u); fill_halo_regions!(v)

        diagnose_w!(w, u, v, grid)
        @test all(interior(w)[:, :, 1] .== 0.0)
    end

    # Determinism: repeated calls produce identical w.
    @testset "Repeated calls produce the same w" begin
        z_faces = [-600.0, -300.0, -100.0, 0.0]
        grid = _make_continuity_grid(4, 4, z_faces)
        u  = Field{Face,   Center, Center}(grid)
        v  = Field{Center, Face,   Center}(grid)
        w1 = Field{Center, Center, Face  }(grid)
        w2 = Field{Center, Center, Face  }(grid)
        set!(u, (x, y, z) -> 0.02 * sin(x))
        set!(v, (x, y, z) -> 0.02 * cos(y))
        fill_halo_regions!(u); fill_halo_regions!(v)

        diagnose_w!(w1, u, v, grid)
        diagnose_w!(w2, u, v, grid)
        @test interior(w1) ≈ interior(w2)
    end
end
