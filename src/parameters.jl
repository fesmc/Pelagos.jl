# Physical constants and tuning parameters for Pelagos.jl / GOLDSTEIN ocean model.
# All values must match the CLIMBER-X Fortran namelist defaults unless noted.
# Units are given in comments.

module Parameters

# ── Reference values ──────────────────────────────────────────────────────────
# rho0 = 1000 matches ocn_par.nml; GOLDSTEIN uses 1000 as Boussinesq reference.
const RHO_0       = 1000.0          # kg m⁻³,  Boussinesq reference density
const CP_OCN      = 3994.0          # J kg⁻¹ K⁻¹, specific heat of seawater
const G           = 9.81            # m s⁻²,   gravitational acceleration
const OMEGA       = 7.2921e-5       # rad s⁻¹, Earth rotation rate
const R_EARTH     = 6.371e6         # m,        Earth radius

# ── Baroclinic velocity solve ──────────────────────────────────────────────────
const R_BC        = 4.0 / 86400.0   # s⁻¹,  baroclinic friction coefficient (4 d⁻¹)
const F_MIN       = 5e-6            # s⁻¹,  Equatorial Coriolis floor (applied to |f|)

# ── Barotropic streamfunction solve ───────────────────────────────────────────
const R_BT_FAC    = 3.0             # –,    near-boundary barotropic friction multiplier
# Base barotropic friction (s⁻¹) — spatially variable; R_BT_FAC applied at boundaries
const R_BT_BASE   = 3.0 / 86400.0  # s⁻¹,  base barotropic friction coefficient

# ── Tracer diffusion ───────────────────────────────────────────────────────────
const K_ISO       = 1500.0          # m² s⁻¹, isopycnal (Redi) diffusivity
const K_GM        = 1500.0          # m² s⁻¹, Gent-McWilliams coefficient
const SLOPE_MAX   = 1e-3            # –,       Gerdes et al. (1991) slope taper limit
const K_CONV      = 100.0           # m² s⁻¹, convective adjustment diffusivity

# Bryan & Lewis (1979) diapycnal diffusivity profile (i_diff_dia=3 in ocn_par.nml)
# κ(z) = K_BG + (K_DEEP-K_BG)·(2/π)·arctan(|z| / K_DEEP_ZREF)
# κ(0) = K_BG (surface minimum); κ(K_DEEP_ZREF) = midpoint; κ(∞) = K_DEEP.
# This form ensures κ ∈ [K_BG, K_DEEP] everywhere.
# Note: exact formula in CLIMBER-X ocn_diff.f90 needs verification against Fortran.
const K_BG       = 1e-5           # m² s⁻¹, surface diapycnal κ (diff_dia_min = 1e-5 m²/s)
const K_DEEP     = 1.5e-4        # m² s⁻¹, deep asymptotic κ   (diff_dia_max = 1.5e-4 m²/s)
const K_DEEP_ZREF = 1000.0       # m,       half-saturation depth (diff_dia_zref = 1000 m)

# ── Geothermal forcing ─────────────────────────────────────────────────────────
const Q_GEO_DEFAULT = 0.05          # W m⁻²,  default uniform geothermal heat flux

# ── Rigid-lid / virtual salinity flux ─────────────────────────────────────────
# (No free parameters; behaviour is determined by the rigid-lid architecture.)

end # module Parameters
