# Frictional-geostrophic baroclinic velocity solver on an Oceananigans C-grid.
#
# Solves at each grid point:
#   -f·v + r_bc·u = Fx     →   u = (r_bc·Fx + f·Fy) / (r_bc² + f²)
#    f·u + r_bc·v = Fy         v = (r_bc·Fy - f·Fx) / (r_bc² + f²)
#
# Grid staggering (C-grid, consistent with Oceananigans):
#   p on T-points  (Center, Center, Center)
#   u on U-points  (Face,   Center, Center)  ← east face of each T-cell
#   v on V-points  (Center, Face,   Center)  ← north face of each T-cell
#
# Pressure gradients:
#   Fx at U-point: exact centred difference (p[i+1,j] - p[i,j]) / Δx
#   Fy at U-point: 4-point average of V-face ∂p/∂y to the U-point
#   Fy at V-point: exact centred difference (p[i,j+1] - p[i,j]) / Δy
#   Fx at V-point: 4-point average of U-face ∂p/∂x to the V-point
#
# The 4-point averaging of the Coriolis cross-term is the standard C-grid
# treatment (see Arakawa & Lamb 1977; Oceananigans momentum equations).
# This is equivalent to the B-grid solve used in GOLDSTEIN in the limit of
# uniform grid spacing, and gives a depth-integrated velocity field that is
# nearly divergence-free on the Oceananigans C-grid.
#
# Reference: Appendix B of Willeit et al. (2022); ocn_baroclinic.f90.

module Baroclinic

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: AbstractGrid, Face, Center
using Oceananigans.Fields: Field
using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ, Δzᵃᵃᶜ
using Oceananigans.Grids: inactive_node, peripheral_node

using ..Parameters: R_BC, F_MIN, RHO_0, OMEGA

export solve_baroclinic!, coriolis_parameter

"""
    coriolis_parameter(lat_deg)

f = 2Ω sin(φ) with the equatorial floor |f| ≥ F_MIN applied.
"""
@inline function coriolis_parameter(lat_deg::Float64)::Float64
    f = 2.0 * OMEGA * sind(lat_deg)
    s = f >= 0.0 ? 1.0 : -1.0
    return s * max(abs(f), F_MIN)
end

function coriolis_parameter(lat_deg::AbstractVector{Float64})
    return [coriolis_parameter(φ) for φ in lat_deg]
end

# ── Pressure gradient helpers ──────────────────────────────────────────────────
# Return Pa m⁻¹ (p is in bar so ×1e5 converts bar→Pa).

@inline function _Fx_at_u(i, j, k, grid, p)
    return ∂xᶠᶜᶜ(i, j, k, grid, p) * 1e5   # at U-face (Face, Center, Center)
end

@inline function _Fy_at_v(i, j, k, grid, p)
    return ∂yᶜᶠᶜ(i, j, k, grid, p) * 1e5   # at V-face (Center, Face, Center)
end

# 4-point average of Fy from surrounding V-faces to the U-face at (i−½, j).
# Oceananigans convention:
#   u at (Face, Center, Center) — u[i,j,k] sits at face (i−½, j),
#                                  between T-cells i−1 and i
#   ∂yᶜᶠᶜ(i,j,k) = (p[i,j,k] − p[i,j−1,k])/Δy — at face (i, j−½),
#                  the south face of T-cell j
# V-faces surrounding U-face (i−½, j):
#   (i−1, j−½) = ∂yᶜᶠᶜ(i−1, j,   k)
#   (i−1, j+½) = ∂yᶜᶠᶜ(i−1, j+1, k)
#   (i,   j−½) = ∂yᶜᶠᶜ(i,   j,   k)
#   (i,   j+½) = ∂yᶜᶠᶜ(i,   j+1, k)
@inline function _Fy_at_u(i, j, k, grid, p)
    return 0.25 * (∂yᶜᶠᶜ(i-1, j,   k, grid, p) +
                   ∂yᶜᶠᶜ(i-1, j+1, k, grid, p) +
                   ∂yᶜᶠᶜ(i,   j,   k, grid, p) +
                   ∂yᶜᶠᶜ(i,   j+1, k, grid, p)) * 1e5
end

# 4-point average of Fx from surrounding U-faces to the V-face at (i, j−½).
# v[i,j,k] sits at face (i, j−½), between T-cells j−1 and j.
# U-faces surrounding V-face (i, j−½):
#   (i−½, j−1) = ∂xᶠᶜᶜ(i,   j−1, k)
#   (i+½, j−1) = ∂xᶠᶜᶜ(i+1, j−1, k)
#   (i−½, j)   = ∂xᶠᶜᶜ(i,   j,   k)
#   (i+½, j)   = ∂xᶠᶜᶜ(i+1, j,   k)
@inline function _Fx_at_v(i, j, k, grid, p)
    return 0.25 * (∂xᶠᶜᶜ(i,   j-1, k, grid, p) +
                   ∂xᶠᶜᶜ(i+1, j-1, k, grid, p) +
                   ∂xᶠᶜᶜ(i,   j,   k, grid, p) +
                   ∂xᶠᶜᶜ(i+1, j,   k, grid, p)) * 1e5
end

# ── Main solver ────────────────────────────────────────────────────────────────

"""
    solve_baroclinic!(u, v, p, tau_x, tau_y, grid)

Compute baroclinic horizontal velocities from the frictional-geostrophic balance.

# Arguments
- `u`    : `Field{Face,   Center, Center}` — zonal velocity output
- `v`    : `Field{Center, Face,   Center}` — meridional velocity output
- `p`    : `Field{Center, Center, Center}` — hydrostatic pressure, bar (halos filled)
- `tau_x`: `Field{Face,   Center, Nothing}` or Matrix — zonal wind stress, N m⁻² (at U-face)
- `tau_y`: `Field{Center, Face,   Nothing}` or Matrix — meridional wind stress (at V-face)
- `grid` : the shared ImmersedBoundaryGrid
"""
function solve_baroclinic!(u    ::Field{Face,   Center, Center},
                           v    ::Field{Center, Face,   Center},
                           p    ::Field{Center, Center, Center},
                           tau_x::AbstractMatrix{Float64},
                           tau_y::AbstractMatrix{Float64},
                           grid ::AbstractGrid)
    Nz = grid.Nz
    Nx = grid.Nx
    Ny = grid.Ny

    u_d = interior(u)
    v_d = interior(v)

    # Top-layer thickness for wind stress application
    dz_top = Δzᵃᵃᶜ(1, 1, Nz, grid)

    rbc  = R_BC
    rho0 = RHO_0

    fill!(v_d, 0.0)   # ensures south (j=1) and north (j=Ny+1) V-faces stay 0

    @inbounds for k in 1:Nz

        # ── u at U-face (i−½, j) = (Face, Center, Center) ─────────────────────
        # Face i lies between T-cells i−1 and i; Periodic in x wraps i=1.
        # Latitude of U-face = latitude of T-row j.
        for j in 1:Ny
            φ_u  = grid.φᵃᶜᵃ[j]
            f_u  = coriolis_parameter(φ_u)
            denom_u = rbc^2 + f_u^2

            for i in 1:Nx
                # Mask using peripheral_node: a U-face is "peripheral" if EITHER
                # neighbouring T-cell is inactive (below bathymetry or land).
                # Required because the FG pressure gradient through a step in
                # bathymetry uses pressures from incompatible depths.
                if peripheral_node(i, j, k, grid, Face(), Center(), Center())
                    u_d[i, j, k] = 0.0
                else
                    Fx_u = -_Fx_at_u(i, j, k, grid, p) / rho0
                    Fy_u = -_Fy_at_u(i, j, k, grid, p) / rho0
                    if k == Nz
                        im1  = mod1(i-1, Nx)
                        τx_u = 0.5 * (tau_x[im1, j] + tau_x[i, j])
                        τy_u = 0.5 * (tau_y[im1, j] + tau_y[i, j])
                        Fx_u += τx_u / (rho0 * dz_top)
                        Fy_u += τy_u / (rho0 * dz_top)
                    end
                    u_d[i, j, k] = (rbc * Fx_u + f_u * Fy_u) / denom_u
                end
            end
        end

        # ── v at V-face (i, j−½) = (Center, Face, Center) ────────────────────
        # v has Ny+1 y-faces; j=1 = south domain wall, j=Ny+1 = north wall.
        # Interior V-faces (j=2..Ny) lie between T-rows j−1 and j.
        for j in 2:Ny
            φ_v  = grid.φᵃᶠᵃ[j]                    # latitude at V-face j
            f_v  = coriolis_parameter(φ_v)
            denom_v = rbc^2 + f_v^2

            for i in 1:Nx
                if peripheral_node(i, j, k, grid, Center(), Face(), Center())
                    v_d[i, j, k] = 0.0
                else
                    Fy_v = -_Fy_at_v(i, j, k, grid, p) / rho0
                    Fx_v = -_Fx_at_v(i, j, k, grid, p) / rho0
                    if k == Nz
                        τx_v = 0.5 * (tau_x[i, j-1] + tau_x[i, j])
                        τy_v = 0.5 * (tau_y[i, j-1] + tau_y[i, j])
                        Fy_v += τy_v / (rho0 * dz_top)
                        Fx_v += τx_v / (rho0 * dz_top)
                    end
                    v_d[i, j, k] = (rbc * Fy_v - f_v * Fx_v) / denom_v
                end
            end
        end
    end

    fill_halo_regions!(u)
    fill_halo_regions!(v)
    return u, v
end

end # module Baroclinic
