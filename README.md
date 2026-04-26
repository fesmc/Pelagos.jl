# Pelagos.jl

[![Build Status](https://github.com/fesmc/Pelagos.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fesmc/Pelagos.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://fesmc.github.io/Pelagos.jl/)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://fesmc.github.io/Pelagos.jl/dev/)

A Julia reimplementation of the **GOLDSTEIN** frictional-geostrophic ocean model,
targeting functional equivalence with the GOLDSTEIN component of
[CLIMBER-X v1.0](https://doi.org/10.5194/gmd-15-5905-2022)
(Willeit et al., GMD 2022).

## What is this?

GOLDSTEIN (Global Ocean Linear Drag Salt and Temperature Equation Integrator) is
a coarse-resolution ocean model used in Earth System Models. Unlike full
Navier-Stokes models, GOLDSTEIN **diagnoses** horizontal velocities algebraically
from a frictional-geostrophic momentum balance — making it computationally
inexpensive while reproducing the large-scale ocean circulation.

**Pelagos.jl** reimplements GOLDSTEIN in idiomatic Julia, using
[Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) as the tracer
transport engine while keeping the velocity solver entirely custom.

## Key features

- All model state lives in **Oceananigans `Field`s on a shared `ImmersedBoundaryGrid`** —
  pressure, velocity, T, S all use the same grid metrics, halos, and immersed-boundary
  masking. No separate index conventions or GridParams; one source of geometry.
- UNESCO (Millero & Poisson 1981) equation of state — not linearised, not TEOS-10
- Frictional-geostrophic baroclinic velocity solve on a C-grid using Oceananigans'
  derivative operators (`∂xᶠᶜᶜ`, `∂yᶜᶠᶜ`) with proper 4-point Coriolis averaging
- Barotropic streamfunction from a 2D elliptic vorticity equation with island
  constraints; corner-ψ formulation makes the depth-integrated transport
  **discretely divergence-free to machine precision**, so `w_top` ≈ 1×10⁻¹⁸ m s⁻¹
- Vertical velocity diagnosed from ∇·**u** = 0 with a custom flux-form continuity
  that respects coastal stepped bathymetry without over-masking
- T/S tracer transport via Oceananigans `PrescribedVelocityFields` with GM/Redi,
  Bryan-Lewis diapycnal mixing, and convective adjustment
- Virtual salinity flux for freshwater forcing (rigid-lid consistent)
- Geothermal heat flux at the ocean floor

## Architecture

```
┌─────────────────────────────────────────────────────────┐
│                     Pelagos.jl                          │
│                                                         │
│  ┌──────────────────┐   u,v,w    ┌─────────────────┐   │
│  │  VelocitySolver  │ ─────────► │  Oceananigans   │   │
│  │  (custom Julia)  │            │  Prescribed     │   │
│  │                  │            │  VelocityFields │   │
│  │  1. UNESCO EOS   │            │  + T,S model    │   │
│  │  2. ∂p/∂z hydro  │            │                 │   │
│  │  3. Baroclinic   │            │  Owns: T, S     │   │
│  │     FG solve     │            │  GM/Redi        │   │
│  │  4. Barotropic   │            │  Bryan-Lewis κ  │   │
│  │     ψ solve      │            │  Convection     │   │
│  │  5. w = ∇·u = 0  │            │                 │   │
│  └──────────────────┘            └─────────────────┘   │
└─────────────────────────────────────────────────────────┘
```

Oceananigans solves no momentum equations and maintains no free surface — both
are handled by Pelagos. It is used as the tracer engine and provides the grid,
`Field`s, halo management, immersed-boundary masking, and finite-difference
operators that the velocity solver builds on.

## Quick start

```julia
using Pelagos

# Build a coupled velocity + tracer ocean model from a CLIMBER-X restart
m = build_ocean_model("path/to/climber-x/restart/pi_cc_open")

# Step it forward (zero wind stress here; pass real τx, τy in production)
tau_x = zeros(72, 36); tau_y = zeros(72, 36)
for step in 1:365
    step_ocean!(m, tau_x, tau_y, 86400.0)   # 1-day step
end
```

Each `step_ocean!` runs the full pipeline: hydrostatic pressure → frictional-
geostrophic baroclinic velocity → barotropic streamfunction + correction →
diagnosed `w` → Oceananigans T/S advance.

## Driver scripts

The `scripts/` directory contains a self-contained 1-year run + plotting
workflow:

```bash
# 1) Run 12 monthly snapshots, T, S, u, v, w, ψ_bt, ρ, p → NetCDF
PELAGOS_NO_WIND=1 julia --project=. scripts/run_1yr.jl

# 2) Make a set of CairoMakie figures from the NetCDF
julia --project=scripts scripts/plot_run.jl scripts/output/pelagos_1yr_monthly.nc
```

Useful environment variables for the run script:
`PELAGOS_RESTART`, `PELAGOS_FORCING`, `PELAGOS_OUTPUT`, `PELAGOS_DT_DAYS`,
`PELAGOS_NO_WIND`. The plotting script auto-instantiates its own
`scripts/Project.toml` (CairoMakie + NCDatasets) on first run.

## Development phases

Development follows a strict sequential order validated against CLIMBER-X Fortran
reference output:

1. ✅ **EOS and pressure** — UNESCO EOS + hydrostatic integration
2. ✅ **Baroclinic velocity** — frictional-geostrophic algebraic solve
3. ✅ **Barotropic solver** — 2D elliptic ψ with island constraints
4. ✅ **Continuity** — diagnose w from ∇·u = 0
5. ✅ **Oceananigans T/S model** — coupled `OceanModel` runs stably for ≥1 year;
   `w_top` at machine precision; tracer drift consistent with diffusion+advection
6. 🔲 **Integration validation** — 100-year spinup vs. CLIMBER-X fixtures

## Running the tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Fixture-based integration tests (Phase 6) require reference NetCDF output from
CLIMBER-X. See `test/fixtures/README.md` for details.

## Reference

Willeit, M., Ganopolski, A., Robinson, A., and Edwards, N. R. (2022).
The Earth system model CLIMBER-X v1.0 – Part 1: Climate model description and
validation. *Geosci. Model Dev.*, 15, 5905–5948.
https://doi.org/10.5194/gmd-15-5905-2022
