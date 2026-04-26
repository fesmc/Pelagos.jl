# Pelagos.jl

**Pelagos.jl** is a Julia reimplementation of the GOLDSTEIN frictional-geostrophic
ocean model, targeting functional equivalence with the GOLDSTEIN component of
CLIMBER-X v1.0 ([Willeit et al., GMD 2022](https://doi.org/10.5194/gmd-15-5905-2022)).

## What is GOLDSTEIN?

GOLDSTEIN (Global Ocean Linear Drag Salt and Temperature Equation Integrator) is a
coarse-resolution (typically 5°×5°) ocean model used within Earth System Models. It is
designed to be computationally inexpensive while reproducing the large-scale thermohaline
circulation and tracer distributions.

Unlike Navier-Stokes models, GOLDSTEIN **diagnoses** horizontal velocities algebraically
from a frictional-geostrophic momentum balance. This makes it far cheaper than
prognostic momentum solvers but limits its application to coarse resolutions where
geostrophic balance is a reasonable assumption.

## Key features

- All state in [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) `Field`s
  on a shared `ImmersedBoundaryGrid` — pressure, velocity, T, S share the same
  grid, halos, and immersed-boundary masking
- UNESCO (Millero & Poisson 1981) equation of state — **not** linearised, **not** TEOS-10
- Frictional-geostrophic baroclinic velocity solve on a C-grid using Oceananigans'
  derivative operators with proper 4-point Coriolis averaging
- Barotropic streamfunction from an elliptic 2D vorticity equation with island
  constraints; corner-ψ formulation makes the depth-integrated transport
  **discretely divergence-free to machine precision**
- Vertical velocity diagnosed from ∇·**u** = 0 with a custom flux-form continuity
  that handles stepped coastal bathymetry correctly
- Tracer (T, S) transport via Oceananigans `PrescribedVelocityFields` with
  GM/Redi isopycnal diffusion and Bryan-Lewis diapycnal mixing
- Virtual salinity flux for freshwater forcing (rigid-lid consistent)
- Geothermal heat flux at the ocean floor

## Quick start

The high-level API is [`build_ocean_model`](@ref) and [`step_ocean!`](@ref):

```julia
using Pelagos

# Initialise a coupled velocity + tracer model from a CLIMBER-X restart
m = build_ocean_model("path/to/climber-x/restart/pi_cc_open")

# Step it forward (1-day timestep)
tau_x = zeros(72, 36); tau_y = zeros(72, 36)
step_ocean!(m, tau_x, tau_y, 86400.0)
```

Each `step_ocean!` runs the full pipeline:

1. Compute hydrostatic pressure from T, S
2. Solve the frictional-geostrophic baroclinic balance for u, v
3. Solve the barotropic streamfunction and apply the C-grid divergence-free correction
4. Diagnose w from ∇·**u** = 0
5. Advance T, S one step in Oceananigans

Lower-level pieces are available too:

```julia
# Compute seawater density at T=10°C, S=35 psu, p=50 bar
ρ = seawater_density(10.0, 35.0, 50.0)

# Coriolis parameter array
f = coriolis_parameter(collect(-80.0:5.0:80.0))
```

## End-to-end example: 1-year monthly run + figures

Two ready-to-run scripts in `scripts/` cover the full workflow:

* `scripts/run_1yr.jl` — drives `OceanModel` for 12 × 30-day months from the
  CLIMBER-X PI restart (or zeroes wind stress with `PELAGOS_NO_WIND=1`),
  writing T, S, u, v, w, ψ\_bt, ρ, p plus the static ocean mask and
  bathymetry to a single NetCDF file.
* `scripts/plot_run.jl` — reads that NetCDF and produces six PNG figures
  (summary, surface maps, zonal-mean depth–latitude sections, zonal-mean
  profiles, MOC, monthly time series) using CairoMakie.

```bash
# Run the diffusion-only baseline (≈ 30 s)
PELAGOS_NO_WIND=1 julia --project=. scripts/run_1yr.jl

# Make figures (CairoMakie installs into scripts/Project.toml on first call)
julia --project=scripts scripts/plot_run.jl scripts/output/pelagos_1yr_monthly.nc
```

## Installation

```julia
using Pkg
Pkg.add("Pelagos")
```

Or from source:

```julia
Pkg.develop(url = "https://github.com/YOUR_ORG/Pelagos.jl")
```

## Citation

If you use Pelagos.jl in published research, please cite the CLIMBER-X paper:

> Willeit, M., Ganopolski, A., Robinson, A., and Edwards, N. R. (2022).
> The Earth system model CLIMBER-X v1.0 – Part 1: Climate model description and
> validation. *Geosci. Model Dev.*, 15, 5905–5948.
> https://doi.org/10.5194/gmd-15-5905-2022
