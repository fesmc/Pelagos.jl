# Hydrostatic pressure integration for GOLDSTEIN ocean model.
#
# Pressure is computed at cell centres by integrating ρ·g downward from the surface.
# The top-surface boundary condition is p_surface = 0 (rigid lid, no sea surface pressure).
# Convention: z is negative downward, cell centres at z[k] < 0, k=1 at surface.
#
# Staggering: pressure at cell centres (T-points), consistent with Fortran ocn_pressure.f90.

module Pressure

using ..UNESCO: seawater_density

export compute_pressure!, compute_pressure

const G   = 9.81    # m s⁻²
const RHO_0 = 1025.0 # kg m⁻³

"""
    compute_pressure!(p, T, S, z_faces, ocean_mask)

Integrate hydrostatic pressure downward from the surface.

# Arguments
- `p`         : output pressure array (nlon, nlat, nz), bar
- `T`         : potential temperature (nlon, nlat, nz), °C
- `S`         : salinity (nlon, nlat, nz), psu
- `z_faces`   : vertical face positions (nz+1,), m, negative downward, z_faces[1]=0
- `ocean_mask`: Bool array (nlon, nlat), true = ocean cell
"""
function compute_pressure!(p        ::AbstractArray{Float64,3},
                           T        ::AbstractArray{Float64,3},
                           S        ::AbstractArray{Float64,3},
                           z_faces  ::AbstractVector{Float64},
                           ocean_mask::AbstractMatrix{Bool})
    nlon, nlat, nz = size(p)

    # Iterate downward from surface (k=1) to bottom (k=nz)
    @inbounds for k in 1:nz
        dz = abs(z_faces[k] - z_faces[k+1])   # layer thickness (positive)
        for j in 1:nlat, i in 1:nlon
            if !ocean_mask[i, j]
                p[i, j, k] = 0.0
                continue
            end
            # Pressure at level k centre = pressure at level k top + ρ·g·(dz/2)
            # For k=1: p_top = 0 (rigid lid)
            p_top = (k == 1) ? 0.0 : p[i, j, k-1] + _half_layer_pressure(T[i,j,k-1], S[i,j,k-1],
                                                                            p[i,j,k-1],
                                                                            z_faces[k-1], z_faces[k])
            # Pressure at cell centre
            rho_k  = seawater_density(T[i,j,k], S[i,j,k], p_top * 0.5)  # first-guess p
            p[i,j,k] = p_top + rho_k * G * (dz * 0.5) / 1e5  # convert Pa → bar
        end
    end
    return p
end

# Pressure increment from level-k centre to the top face of level k+1
@inline function _half_layer_pressure(T::Float64, S::Float64, p::Float64,
                                      z_top::Float64, z_bot::Float64)::Float64
    dz   = abs(z_top - z_bot)
    rho  = seawater_density(T, S, p)
    return rho * G * (dz * 0.5) / 1e5   # bar
end

"""
    compute_pressure(T, S, z_faces, ocean_mask)

Allocating version of `compute_pressure!`. Returns pressure array (nlon, nlat, nz) in bar.
"""
function compute_pressure(T        ::AbstractArray{Float64,3},
                          S        ::AbstractArray{Float64,3},
                          z_faces  ::AbstractVector{Float64},
                          ocean_mask::AbstractMatrix{Bool})
    p = zeros(Float64, size(T))
    compute_pressure!(p, T, S, z_faces, ocean_mask)
    return p
end

end # module Pressure
