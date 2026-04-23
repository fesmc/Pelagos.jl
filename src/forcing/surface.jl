# Surface boundary condition handling for Pelagos.jl.
#
# Implements:
#   - Heat flux as a Neumann condition on T (surface layer flux divergence).
#   - Virtual salinity flux (rigid-lid approximation for freshwater forcing):
#       FS = FW · S_local / H_top
#     where FW is the net freshwater flux (E − P − R) in m s⁻¹.
#   - Global correction: ensure globally integrated virtual salinity flux is zero
#     each timestep (Appendix B4, Willeit et al. 2022).
#   - Wind stress τˣ, τʸ (N m⁻²): passed directly to the baroclinic solver.
#
# None of these functions interact with Oceananigans internals; they compute
# boundary condition arrays that are then applied via Oceananigans FluxBoundaryCondition
# or directly in the velocity solve.

module Surface

export virtual_salinity_flux!, apply_global_salt_correction!

"""
    virtual_salinity_flux!(FS, FW, S_surface, H_top, ocean_mask)

Compute the virtual salinity flux FS (psu m s⁻¹) at the ocean surface.

# Arguments
- `FS`         : output virtual salinity flux (nlon, nlat), psu·m·s⁻¹
- `FW`         : net freshwater flux = E − P − R (nlon, nlat), m s⁻¹ (positive = into ocean)
- `S_surface`  : surface salinity (nlon, nlat), psu
- `H_top`      : top layer thickness (m), scalar
- `ocean_mask` : Bool (nlon, nlat)
"""
function virtual_salinity_flux!(FS        ::AbstractMatrix{Float64},
                                FW        ::AbstractMatrix{Float64},
                                S_surface ::AbstractMatrix{Float64},
                                H_top     ::Float64,
                                ocean_mask::AbstractMatrix{Bool})
    @. FS = ifelse(ocean_mask, FW * S_surface / H_top, 0.0)
    return FS
end

"""
    apply_global_salt_correction!(FS, ocean_mask, area)

Adjust the virtual salinity flux field so its globally integrated value is zero.
This preserves total ocean salinity (Appendix B4 of Willeit et al. 2022).

# Arguments
- `FS`         : virtual salinity flux (nlon, nlat), psu·m·s⁻¹ — modified in place
- `ocean_mask` : Bool (nlon, nlat)
- `area`       : cell area (nlon, nlat), m²
"""
function apply_global_salt_correction!(FS        ::AbstractMatrix{Float64},
                                       ocean_mask::AbstractMatrix{Bool},
                                       area      ::AbstractMatrix{Float64})
    total_flux = 0.0
    total_area = 0.0
    nlon, nlat = size(FS)
    @inbounds for j in 1:nlat, i in 1:nlon
        ocean_mask[i, j] || continue
        total_flux += FS[i, j] * area[i, j]
        total_area += area[i, j]
    end
    total_area > 0.0 || return FS
    correction = total_flux / total_area
    @inbounds for j in 1:nlat, i in 1:nlon
        ocean_mask[i, j] || continue
        FS[i, j] -= correction
    end
    return FS
end

"""
    heat_flux_bc(QH, H_top, rho0, cp)

Convert surface heat flux QH (W m⁻²) to a temperature tendency (K s⁻¹)
for the top layer. Positive QH = into ocean.

Returns a (nlon, nlat) array of ∂T/∂t forcing due to surface heating.
"""
function heat_flux_bc(QH   ::AbstractMatrix{Float64},
                      H_top::Float64,
                      rho0 ::Float64 = 1025.0,
                      cp   ::Float64 = 3994.0)::Matrix{Float64}
    return @. QH / (rho0 * cp * H_top)
end

end # module Surface
