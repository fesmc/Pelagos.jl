# Hydrostatic pressure integration on an Oceananigans grid.
#
# p=0 at the surface (rigid lid). Integrates downward from k=Nz (shallowest)
# to k=1 (deepest) — Oceananigans' bottom-to-top index convention.
# Pressure is in bar (1 bar = 1e5 Pa).
# After integration, fill_halo_regions! is called so p can be used directly
# with Oceananigans derivative operators in the baroclinic solver.

module Pressure

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: AbstractGrid
using Oceananigans.Fields: Field, Center
using Oceananigans.Operators: Δzᵃᵃᶜ
using Oceananigans.Grids: inactive_node

using ..UNESCO: seawater_density

export compute_pressure!

const _G = 9.81       # m s⁻²

"""
    compute_pressure!(p, T, S, grid)

Integrate hydrostatic pressure downward into the water column.

`p`, `T`, `S` are `Field{Center, Center, Center}` on `grid`.
Uses Oceananigans k-convention: k=Nz = surface layer, k=1 = deepest layer.
Immersed (below-bathymetry) cells are set to p = pressure of the nearest
active cell above them, so horizontal pressure gradients at the seafloor
boundary remain small (avoids spurious velocities from the FG solve).
Returns `p` with halos filled.
"""
function compute_pressure!(p    ::Field{Center, Center, Center},
                           T    ::Field{Center, Center, Center},
                           S    ::Field{Center, Center, Center},
                           grid ::AbstractGrid)
    Nz = grid.Nz
    Nx = grid.Nx
    Ny = grid.Ny

    p_d = interior(p)
    T_d = interior(T)
    S_d = interior(S)

    # Surface layer (k = Nz): pressure at the centre of the top cell.
    # p_surface = 0 (rigid lid), so p[Nz] = ρ·g·(Δz/2).
    @inbounds for j in 1:Ny, i in 1:Nx
        if inactive_node(i, j, Nz, grid, Center(), Center(), Center())
            p_d[i, j, Nz] = 0.0
            continue
        end
        Tk  = T_d[i, j, Nz]
        Sk  = max(S_d[i, j, Nz], 0.0)
        dzk = Δzᵃᵃᶜ(i, j, Nz, grid)
        ρ   = seawater_density(Tk, Sk, 0.0)
        p_d[i, j, Nz] = ρ * _G * (dzk * 0.5) / 1e5
    end

    # Integrate downward: k = Nz-1 → 1.
    @inbounds for k in Nz-1:-1:1
        for j in 1:Ny, i in 1:Nx
            if inactive_node(i, j, k, grid, Center(), Center(), Center())
                # For immersed cells copy the pressure from the layer above
                # so horizontal pressure gradients at the seafloor stay small.
                p_d[i, j, k] = p_d[i, j, k+1]
                continue
            end
            dz_above = Δzᵃᵃᶜ(i, j, k+1, grid)
            dz_below = Δzᵃᵃᶜ(i, j, k,   grid)

            Tk_above = T_d[i, j, k+1];  Sk_above = max(S_d[i, j, k+1], 0.0)
            pk_above = p_d[i, j, k+1]
            Tk       = T_d[i, j, k];    Sk       = max(S_d[i, j, k],   0.0)

            ρ_above = seawater_density(Tk_above, Sk_above, pk_above)
            p_face  = pk_above + ρ_above * _G * (dz_above * 0.5) / 1e5

            ρ_below = seawater_density(Tk, Sk, p_face)
            p_d[i, j, k] = p_face + ρ_below * _G * (dz_below * 0.5) / 1e5
        end
    end

    # 2-D land-column ghost fill: copy pressure from an adjacent ocean column
    # so ∂p/∂x and ∂p/∂y at ocean-land boundaries stay small.  Without this,
    # the pressure gradient at a wall gives (p_ocean − 0)/Δx which would
    # produce catastrophic velocities in the FG solve.
    @inbounds for j in 1:Ny, i in 1:Nx
        # Is this column fully land? (top cell is inactive)
        inactive_node(i, j, Nz, grid, Center(), Center(), Center()) || continue
        ip1 = mod1(i+1, Nx); im1 = mod1(i-1, Nx)
        jp1 = min(j+1, Ny);  jm1 = max(j-1, 1)
        src_i, src_j = 0, 0
        for (ci, cj) in ((ip1, j), (im1, j), (i, jp1), (i, jm1))
            if !inactive_node(ci, cj, Nz, grid, Center(), Center(), Center())
                src_i, src_j = ci, cj
                break
            end
        end
        if src_i != 0
            for k in 1:Nz
                p_d[i, j, k] = p_d[src_i, src_j, k]
            end
        end
    end

    fill_halo_regions!(p)
    return p
end

end # module Pressure
