"""
    Pelagos

Julia reimplementation of the GOLDSTEIN frictional-geostrophic ocean model,
targeting functional equivalence with the GOLDSTEIN component of CLIMBER-X v1.0
(Willeit et al., GMD 2022, https://doi.org/10.5194/gmd-15-5905-2022).

Architecture:
- Custom velocity solver (frictional-geostrophic + barotropic ψ) lives in src/velocity/
- Oceananigans is used **only** as the tracer transport engine via PrescribedVelocityFields
- UNESCO (Millero & Poisson 1981) equation of state implemented from scratch in src/eos/

Development proceeds in phases; do not use modules from later phases before
the tests for earlier phases pass.
"""
module Pelagos

# ── Physical constants and parameters ──────────────────────────────────────────
include("parameters.jl")
using .Parameters

# ── Equation of state and hydrostatic pressure ─────────────────────────────────
include("eos/UNESCO.jl")
include("eos/pressure.jl")
using .UNESCO
using .Pressure

# ── Velocity solver ────────────────────────────────────────────────────────────
include("velocity/islands.jl")
include("velocity/baroclinic.jl")
include("velocity/barotropic.jl")
include("velocity/continuity.jl")
using .Islands
using .Baroclinic
using .Barotropic
using .Continuity

# ── Grid utilities ─────────────────────────────────────────────────────────────
include("grid/bathymetry.jl")
include("grid/masks.jl")
include("grid/setup.jl")
using .Bathymetry
using .Masks
using .GridSetup

# ── Forcing ────────────────────────────────────────────────────────────────────
include("forcing/surface.jl")
include("forcing/geothermal.jl")
using .Surface
using .Geothermal

# ── Tracer model (Oceananigans) ────────────────────────────────────────────────
include("tracers/diffusion.jl")
include("tracers/convection.jl")
include("tracers/setup.jl")
using .Diffusion
using .Convection
using .TracerSetup

# ── I/O ────────────────────────────────────────────────────────────────────────
include("io/output.jl")
include("io/restarts.jl")
using .Output
using .Restarts

# ── Public API re-exports ──────────────────────────────────────────────────────
# EOS
export seawater_density, seawater_density_insitu, seawater_density!
# Pressure
export compute_pressure, compute_pressure!
# Velocity
export solve_baroclinic!, coriolis_parameter
export build_barotropic_solver, solve_barotropic!, compute_rbt
export diagnose_w, diagnose_w!
# Islands
export detect_islands, IslandInfo
# Grid
export smooth_bathymetry, load_bathymetry_nc
export build_ocean_mask, build_velocity_mask, build_tracer_mask
export build_climberx_grid
# Forcing
export virtual_salinity_flux!, apply_global_salt_correction!, heat_flux_bc
export geothermal_tendency, uniform_geothermal
# Tracers
export bryan_lewis_kappa, gm_redi_closure, diapycnal_closure
export convective_adjustment_closure
export build_tracer_model, update_velocities!
# I/O
export write_snapshot, write_restart, read_restart
export load_climberx_restart, load_climberx_forcing, bgrid_u_to_cgrid, bgrid_v_to_cgrid

end # module Pelagos
