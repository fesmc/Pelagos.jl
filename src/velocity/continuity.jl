# Vertical velocity diagnosis from the continuity equation (∇·u = 0).
#
# After u and v are computed from the frictional-geostrophic solve, w is
# determined by integrating the continuity equation upward from the ocean floor
# (w = 0 at bathymetry):
#
#   w[k] = w[k+1] - (∂u/∂x + ∂v/∂y)·Δz[k]    (k from bottom to surface)
#
# The rigid-lid constraint requires w[1] (top of surface layer) ≈ 0.
# If |w_surface| > 1e-12 m/s, there is a barotropic solve error.
#
# Grid convention: w on W-points (vertical cell faces), positive upward.
# Metric factors on the sphere: ∂u/∂x → (1/a·cosφ)·∂u/∂λ, ∂v/∂y → (1/a)·∂v/∂φ.

module Continuity

using ..Parameters: R_EARTH

export diagnose_w!, diagnose_w

"""
    diagnose_w!(w, u, v, dx, dy_v, cos_lat, dz, ocean_mask; debug=false)

Compute vertical velocity from depth-integrated continuity.

# Arguments
- `w`         : output vertical velocity (nlon, nlat, nz+1), m s⁻¹, on W-points
                 w[:,:,1] = surface (rigid lid → ≈0), w[:,:,nz+1] = bottom = 0
- `u`         : zonal velocity (nlon, nlat, nz), m s⁻¹
- `v`         : meridional velocity (nlon, nlat, nz), m s⁻¹
- `dx`        : zonal grid spacing (nlat,), m (at u-point latitudes)
- `dy_v`      : meridional grid spacing at v-points, m (scalar, uniform assumed)
- `cos_lat`   : cos(latitude) at T-point centres (nlat,)
- `dz`        : layer thicknesses (nz,), m
- `ocean_mask`: Bool (nlon, nlat)
- `debug`     : if true, assert that surface w is below 1e-10 m s⁻¹
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
    fill!(w, 0.0)

    # Integrate from bottom (k = nz) upward to surface (k = 1)
    # w[:,:,nz+1] = 0  (bottom, already set)
    @inbounds for k in nz:-1:1
        dzk  = dz[k]
        for j in 1:nlat
            dxj    = dx[j]
            cosφ_j = cos_lat[j]
            # cosφ at v-face (between j and j+1): approximate as average
            cosφ_N = j < nlat ? 0.5*(cos_lat[j] + cos_lat[j+1]) : 0.0
            cosφ_S = j >    1 ? 0.5*(cos_lat[j] + cos_lat[j-1]) : 0.0

            for i in 1:nlon
                ocean_mask[i, j] || continue
                ip1 = mod1(i+1, nlon)

                # ∂u/∂x: (u[i+½,j,k] - u[i-½,j,k]) / dx
                # u is stored on east face: u[i,j,k] ≡ u at east face of cell (i,j)
                du_dx = (u[i, j, k] - u[mod1(i-1,nlon), j, k]) / dxj

                # ∂(v·cosφ)/∂y: metric factor for sphere
                vN = j < nlat ? v[i, j, k]   : 0.0
                vS = j >    1 ? v[i, j-1, k] : 0.0
                # v[i,j,k] is on the north face of cell (i,j)
                dv_dy = (vN * cosφ_N - vS * cosφ_S) / (dy_v * cosφ_j)

                # Upward integration: w(top of k) = w(bottom of k) + divergence·Δz
                w[i, j, k] = w[i, j, k+1] + (du_dx + dv_dy) * dzk
            end
        end
    end

    if debug
        for j in 1:nlat, i in 1:nlon
            ocean_mask[i, j] || continue
            @assert abs(w[i, j, 1]) < 1e-10 "Rigid-lid violation: w_surface[$(i),$(j)] = $(w[i,j,1])"
        end
    end
    return w
end

"""
    diagnose_w(u, v, dx, dy_v, cos_lat, dz, ocean_mask; debug=false)

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
