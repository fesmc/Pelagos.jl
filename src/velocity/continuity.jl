# Vertical velocity diagnosis from the continuity equation (∇·u = 0).
#
# After u and v are computed from the frictional-geostrophic solve, w is
# determined by integrating the continuity equation upward from the ocean floor
# (w = 0 at bathymetry):
#
#   w[k+1] = w[k] - (∂u/∂x + ∂v/∂y)·Δz[k]   (k = 1 → nz, bottom-to-top)
#
# Array conventions (same as Oceananigans and the rest of this project):
#   - Layer k=1 is the deepest (bottom), k=nz is the shallowest (surface).
#   - w is on W-points (vertical cell faces), positive upward.
#   - w[:,:,1]    = 0   at bathymetry  (bottom face of layer 1)
#   - w[:,:,nz+1] ≈ 0   at sea surface (rigid-lid constraint)
#
# Metric factors on the sphere:
#   ∂u/∂x → (u_east - u_west) / (a·cosφ·Δλ)
#   ∂(v·cosφ)/∂y / cosφ → (v_N·cosφ_N - v_S·cosφ_S) / (a·Δφ·cosφ)

module Continuity

using ..Parameters: R_EARTH

export diagnose_w!, diagnose_w

"""
    diagnose_w!(w, u, v, dx, dy_v, cos_lat, dz, ocean_mask; debug=false)

Compute vertical velocity from depth-integrated continuity.

# Arguments
- `w`         : output (nlon, nlat, nz+1), m s⁻¹, on W-points
                 w[:,:,1] = 0 (bathymetry); w[:,:,nz+1] ≈ 0 (rigid lid)
- `u`         : zonal velocity (nlon, nlat, nz), m s⁻¹, east-face C-grid, k=1=deepest
- `v`         : meridional velocity (nlon, nlat, nz), m s⁻¹, north-face of each T-cell,
                 k=1=deepest; v[:,j,:] is the velocity at the north face of T-cell j
- `dx`        : zonal grid spacing (nlat,), m, at T-cell centre latitudes
- `dy_v`      : meridional grid spacing, m (scalar, assumed uniform)
- `cos_lat`   : cos(latitude) at T-cell centres (nlat,)
- `dz`        : layer thicknesses (nz,), m, dz[1]=deepest layer thickness
- `ocean_mask`: Bool (nlon, nlat)
- `debug`     : if true, assert surface w < 1×10⁻¹⁰ m s⁻¹ (rigid-lid check)
"""
function diagnose_w!(w         ::AbstractArray{Float64,3},
                     u         ::AbstractArray{Float64,3},
                     v         ::AbstractArray{Float64,3},
                     dx        ::AbstractVector{Float64},
                     dy_v      ::Float64,
                     cos_lat   ::AbstractVector{Float64},
                     dz        ::AbstractVector{Float64},
                     ocean_mask::AbstractMatrix{Bool};
                     debug::Bool = false)
    nlon, nlat, nz = size(u)
    fill!(w, 0.0)   # w[:,:,1] = 0 at bathymetry; all faces zeroed for land

    # Integrate upward from k=1 (deepest layer) to k=nz (shallowest layer).
    # Continuity: ∂u/∂x + ∂v/∂y + ∂w/∂z = 0
    #   → w_top - w_bottom = -(du_dx + dv_dy) · Δz
    #   → w[i,j,k+1] = w[i,j,k] - (du_dx + dv_dy) · dz[k]
    @inbounds for k in 1:nz
        dzk = dz[k]
        for j in 1:nlat
            dxj    = dx[j]
            cosφ_j = cos_lat[j]
            cosφ_N = j < nlat ? 0.5*(cos_lat[j] + cos_lat[j+1]) : 0.0
            cosφ_S = j >    1 ? 0.5*(cos_lat[j] + cos_lat[j-1]) : 0.0

            for i in 1:nlon
                ocean_mask[i, j] || continue
                # u[i,j,k] is the east-face velocity of cell (i,j).
                # West face of cell i = east face of cell i-1 (periodic).
                du_dx = (u[i, j, k] - u[mod1(i-1, nlon), j, k]) / dxj

                # v[i,j,k] is the north-face velocity of T-cell j.
                # South face of cell j = north face of cell j-1 = v[i,j-1,k].
                vN = j < nlat ? v[i, j,   k] : 0.0
                vS = j >    1 ? v[i, j-1, k] : 0.0
                dv_dy = (vN * cosφ_N - vS * cosφ_S) / (dy_v * cosφ_j)

                w[i, j, k+1] = w[i, j, k] - (du_dx + dv_dy) * dzk
            end
        end
    end

    if debug
        for j in 1:nlat, i in 1:nlon
            ocean_mask[i, j] || continue
            @assert abs(w[i, j, nz+1]) < 1e-10 "Rigid-lid violation: w_surface[$i,$j] = $(w[i,j,nz+1])"
        end
    end
    return w
end

"""
    diagnose_w(u, v, dx, dy_v, cos_lat, dz, ocean_mask; debug=false) -> w

Allocating version of `diagnose_w!`. Returns w (nlon, nlat, nz+1).
"""
function diagnose_w(u         ::AbstractArray{Float64,3},
                    v         ::AbstractArray{Float64,3},
                    dx        ::AbstractVector{Float64},
                    dy_v      ::Float64,
                    cos_lat   ::AbstractVector{Float64},
                    dz        ::AbstractVector{Float64},
                    ocean_mask::AbstractMatrix{Bool};
                    debug::Bool = false)
    nlon, nlat, nz = size(u)
    w = zeros(Float64, nlon, nlat, nz+1)
    diagnose_w!(w, u, v, dx, dy_v, cos_lat, dz, ocean_mask; debug)
    return w
end

end # module Continuity
