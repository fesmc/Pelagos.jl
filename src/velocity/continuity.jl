# Vertical velocity from horizontal continuity: ∇·u = 0.
#
# Integrates upward from the seafloor (w=0 at bathymetry):
#   w[k+1] = w[k] − div_h(u,v) · Δz[k]
#
# Custom flux-form divergence at T-cell (i, j, k):
#     div_h = (1/V) · [Ax(i+1) u(i+1) − Ax(i) u(i) + Ay(j+1) v(j+1) − Ay(j) v(j)]
#
# We deliberately DO NOT use Oceananigans' div_xyᶜᶜᶜ on the ImmersedBoundaryGrid
# because that operator applies a `conditional_δ` that zeroes the *entire*
# divergence at a T-cell as soon as any neighbouring face is immersed.  For
# rigid-lid continuity on a coastal column the open-face fluxes still need to
# contribute; the closed-face contribution is already zero because the
# baroclinic/barotropic solvers set u, v = 0 on closed faces.  Using the
# underlying-grid metrics directly avoids the conditional zeroing while still
# giving correct fluxes (Ax · 0 = 0 at closed faces).
#
# Grid conventions (Oceananigans):
#   k=1 = deepest layer, k=Nz = shallowest layer
#   w on (Center, Center, Face); w[:,:,1] = 0 at bathymetry, w[:,:,Nz+1] ≈ 0
#   u on (Face, Center, Center), v on (Center, Face, Center) — C-grid faces

module Continuity

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: AbstractGrid, Face, Center
using Oceananigans.Fields: Field
using Oceananigans.Operators: Axᶠᶜᶜ, Ayᶜᶠᶜ, Azᶜᶜᶠ, Δzᵃᵃᶜ
using Oceananigans.Grids: inactive_node

export diagnose_w!

# Underlying (un-masked) horizontal flux divergence at T-cell (i, j, k).
# Uses underlying_grid metrics so the conditional-δ on IBG does not interfere.
# Indexes u, v through the Field interface so periodic halos are accessible
# (a raw `interior(u)` array has size (Nx, Ny, Nz) and would error at i=Nx
# when reading u[i+1, j, k] across the dateline).
@inline function _flux_div_h(i, j, k, ug, u, v)
    @inbounds begin
        flux_e = Axᶠᶜᶜ(i+1, j, k, ug) * u[i+1, j, k]
        flux_w = Axᶠᶜᶜ(i,   j, k, ug) * u[i,   j, k]
        flux_n = Ayᶜᶠᶜ(i, j+1, k, ug) * v[i, j+1, k]
        flux_s = Ayᶜᶠᶜ(i, j,   k, ug) * v[i, j,   k]
    end
    return flux_e - flux_w + flux_n - flux_s
end

"""
    diagnose_w!(w, u, v, grid)

Integrate continuity upward from the seafloor to diagnose w.

The surface value `w[:, :, Nz+1]` is the rigid-lid residual; it is ≈ 0 only if
the depth-integrated horizontal divergence of (u, v) is ≈ 0.
"""
function diagnose_w!(w   ::Field{Center, Center, Face},
                     u   ::Field{Face,   Center, Center},
                     v   ::Field{Center, Face,   Center},
                     grid::AbstractGrid)
    Nx = grid.Nx
    Ny = grid.Ny
    Nz = grid.Nz
    ug = grid isa Oceananigans.ImmersedBoundaries.ImmersedBoundaryGrid ? grid.underlying_grid : grid

    w_d = interior(w)

    fill!(w, 0.0)   # w[:,:,1] = 0 at bathymetry

    @inbounds for k in 1:Nz
        for j in 1:Ny, i in 1:Nx
            inactive_node(i, j, k, grid, Center(), Center(), Center()) && continue
            # w[k+1] · A(k+1) − w[k] · A(k) = −flux_div_h    (volume balance)
            # On LLGrid without partial cells, A is uniform with k → divide by A_top.
            A_top = Azᶜᶜᶠ(i, j, k+1, ug)
            w_d[i, j, k+1] = w_d[i, j, k] - _flux_div_h(i, j, k, ug, u, v) / A_top
        end
    end

    fill_halo_regions!(w)
    return w
end

end # module Continuity
