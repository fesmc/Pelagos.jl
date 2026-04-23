# Geothermal heat flux bottom boundary condition for Pelagos.jl.
#
# Following Lucazeau (2019), the geothermal heat flux is applied as a Neumann
# condition on temperature at the ocean floor.  The default is a spatially
# uniform flux Q_GEO = 0.05 W m⁻² (Appendix B of Willeit et al. 2022).
# A spatially varying field from the Lucazeau (2019) dataset can be loaded
# from a NetCDF file when available.

module Geothermal

using ..Parameters: Q_GEO_DEFAULT

export geothermal_tendency

"""
    geothermal_tendency(Q_geo, H_bottom, ocean_mask; rho0=1025.0, cp=3994.0)

Convert geothermal heat flux (W m⁻²) to a temperature tendency (K s⁻¹)
for the bottom ocean layer.

# Arguments
- `Q_geo`     : geothermal heat flux (nlon, nlat), W m⁻²
- `H_bottom`  : bottom layer thickness (m)
- `ocean_mask`: Bool (nlon, nlat)
- `rho0`      : reference density, kg m⁻³
- `cp`        : specific heat, J kg⁻¹ K⁻¹
"""
function geothermal_tendency(Q_geo     ::AbstractMatrix{Float64},
                             H_bottom  ::Float64,
                             ocean_mask::AbstractMatrix{Bool};
                             rho0::Float64 = 1025.0,
                             cp  ::Float64 = 3994.0)::Matrix{Float64}
    dT = similar(Q_geo)
    @. dT = ifelse(ocean_mask, Q_geo / (rho0 * cp * H_bottom), 0.0)
    return dT
end

"""
    uniform_geothermal(nlon, nlat)

Return a uniform geothermal flux field using the default value Q_GEO_DEFAULT.
"""
function uniform_geothermal(nlon::Int, nlat::Int)::Matrix{Float64}
    return fill(Q_GEO_DEFAULT, nlon, nlat)
end

end # module Geothermal
