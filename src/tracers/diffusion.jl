# Tracer diffusion closures for Pelagos.jl.
#
# GM/Redi: delegates to Oceananigans' IsopycnalSkewSymmetricDiffusivity with
#   slope tapering following Gerdes et al. (1991).
#
# Diapycnal: Bryan & Lewis (1979) depth-dependent profile, implemented as
#   a VerticalScalarDiffusivity whose κ depends on z.

module Diffusion

using Oceananigans.TurbulenceClosures
using ..Parameters: K_ISO, K_GM, SLOPE_MAX, K_BG, K_DEEP, K_DEEP_ALPHA, K_DEEP_ZREF

export gm_redi_closure, diapycnal_closure, bryan_lewis_kappa

"""
    bryan_lewis_kappa(z)

Bryan & Lewis (1979) diapycnal diffusivity profile.

κ_d(z) = K_BG + (K_DEEP − K_BG) · (2/π) · arctan(α·(|z| − z_ref))

Returns κ in m² s⁻¹. `z` is depth in metres, negative downward.
"""
@inline function bryan_lewis_kappa(z::Float64)::Float64
    depth = abs(z)
    return K_BG + (K_DEEP - K_BG) * (2.0/π) * atan(K_DEEP_ALPHA * (depth - K_DEEP_ZREF))
end

"""
    gm_redi_closure()

Oceananigans closure implementing the GM + Redi scheme with Gerdes et al. (1991)
slope tapering. κ_iso = κ_gm = 1500 m² s⁻¹, slope_max = 1e-3.
"""
function gm_redi_closure()
    return IsopycnalSkewSymmetricDiffusivity(
        κ_skew     = K_GM,
        κ_symmetric = K_ISO,
        slope_limiter = FluxTapering(SLOPE_MAX),
    )
end

"""
    diapycnal_closure()

Oceananigans vertical diffusivity closure using the Bryan & Lewis (1979) profile.
Returns a `VerticalScalarDiffusivity` with depth-dependent κ.
"""
function diapycnal_closure()
    # Wrap the Bryan-Lewis profile as an Oceananigans diffusivity function.
    # Oceananigans VerticalScalarDiffusivity accepts a function κ(x,y,z,t).
    kappa_func(x, y, z, t) = bryan_lewis_kappa(z)
    return VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), κ = kappa_func)
end

end # module Diffusion
