# Frictional-geostrophic baroclinic velocity solve.
#
# Solves the 2×2 algebraic system at each grid column and level:
#   -f·v + r_bc·u = Fx
#    f·u + r_bc·v = Fy
# → u = (r_bc·Fx + f·Fy) / (r_bc² + f²)
#    v = (r_bc·Fy - f·Fx) / (r_bc² + f²)
#
# where Fx = -(1/ρ₀)·∂p/∂x + τˣ/h  (top layer only)
#       Fy = -(1/ρ₀)·∂p/∂y + τʸ/h
#
# The Equatorial f-floor (F_MIN) is applied consistently here and in barotropic.jl.
# Vectorised over (nlon, nlat) at each level; no explicit horizontal loops in Julia.
#
# Grid convention: pressure on T-points (cell centres), u on U-points (east face),
# v on V-points (north face), consistent with Fortran ocn_baroclinic.f90.

module Baroclinic

using ..Parameters: R_BC, F_MIN, RHO_0, OMEGA, R_EARTH

export solve_baroclinic!, coriolis_parameter

"""
    coriolis_parameter(lat_deg)

Coriolis parameter f = 2Ω·sin(φ) with Equatorial floor applied.
`lat_deg` may be a scalar or array of latitudes in degrees.
"""
@inline function coriolis_parameter(lat_deg::Float64)::Float64
    f    = 2.0 * OMEGA * sind(lat_deg)
    # sign(0.0) == 0 in Julia; treat equator as Northern Hemisphere for the floor.
    s    = f >= 0.0 ? 1.0 : -1.0
    return s * max(abs(f), F_MIN)
end

# Vectorised form over latitude array
function coriolis_parameter(lat_deg::AbstractVector{Float64})
    return [coriolis_parameter(φ) for φ in lat_deg]
end

"""
    pressure_gradient_x(p, i, j, dx)

Zonal pressure gradient at U-point (i+½, j).
p is on T-points, dx is the zonal grid spacing at this latitude (m).
Uses centred difference across the T–T pair straddling the U-point.
"""
@inline function pressure_gradient_x(p::AbstractArray{Float64,3},
                                     i::Int, j::Int, k::Int,
                                     dx::Float64)::Float64
    nlon = size(p, 1)
    ip1  = mod1(i + 1, nlon)   # periodic in longitude
    return (p[ip1, j, k] - p[i, j, k]) / dx
end

"""
    pressure_gradient_y(p, i, j, dy)

Meridional pressure gradient at V-point (i, j+½). No meridional periodicity.
Returns 0 at j = nlat (northern boundary).
"""
@inline function pressure_gradient_y(p::AbstractArray{Float64,3},
                                     i::Int, j::Int, k::Int,
                                     dy::Float64)::Float64
    nlat = size(p, 2)
    j == nlat && return 0.0
    return (p[i, j+1, k] - p[i, j, k]) / dy
end

"""
    solve_baroclinic!(u, v, p, tau_x, tau_y, f, dx, dy, dz, ocean_mask)

Compute baroclinic horizontal velocities from the frictional-geostrophic balance.

# Arguments
- `u`          : output zonal velocity (nlon, nlat, nz), m s⁻¹, on U-points
- `v`          : output meridional velocity (nlon, nlat, nz), m s⁻¹, on V-points
- `p`          : hydrostatic pressure (nlon, nlat, nz), bar, on T-points
- `tau_x`      : zonal wind stress (nlon, nlat), N m⁻², on U-points
- `tau_y`      : meridional wind stress (nlon, nlat), N m⁻², on V-points
- `f`          : Coriolis parameter (nlat,), s⁻¹, with F_MIN floor applied
- `dx`         : zonal grid spacing (nlat,), m, at U-point latitudes
- `dy`         : meridional grid spacing, scalar or (nlat,), m
- `dz`         : layer thicknesses (nz,), m
- `ocean_mask` : Bool (nlon, nlat), true = ocean
"""
function solve_baroclinic!(u        ::AbstractArray{Float64,3},
                           v        ::AbstractArray{Float64,3},
                           p        ::AbstractArray{Float64,3},
                           tau_x    ::AbstractMatrix{Float64},
                           tau_y    ::AbstractMatrix{Float64},
                           f        ::AbstractVector{Float64},
                           dx       ::AbstractVector{Float64},
                           dy       ::Float64,
                           dz       ::AbstractVector{Float64},
                           ocean_mask::AbstractMatrix{Bool})
    nlon, nlat, nz = size(p)
    rbc  = R_BC
    rho0 = RHO_0

    @inbounds for k in 1:nz
        h_top = dz[1]   # surface layer thickness for wind stress (only used at k=1)
        for j in 1:nlat
            fj   = f[j]
            denom = rbc^2 + fj^2
            dxj  = dx[j]

            for i in 1:nlon
                if !ocean_mask[i, j]
                    u[i, j, k] = 0.0
                    v[i, j, k] = 0.0
                    continue
                end

                # Pressure gradient forcing (bar m⁻¹ → Pa m⁻¹ → N m⁻³; ×1e5)
                dpdx = pressure_gradient_x(p, i, j, k, dxj) * 1e5  # Pa m⁻¹
                dpdy = pressure_gradient_y(p, i, j, k, dy)  * 1e5

                Fx = -dpdx / rho0
                Fy = -dpdy / rho0

                # Wind stress forcing: only in top layer
                if k == 1
                    Fx += tau_x[i, j] / (rho0 * h_top)
                    Fy += tau_y[i, j] / (rho0 * h_top)
                end

                u[i, j, k] = (rbc * Fx + fj * Fy) / denom
                v[i, j, k] = (rbc * Fy - fj * Fx) / denom
            end
        end
    end
    return u, v
end

end # module Baroclinic
