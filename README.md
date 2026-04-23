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

- UNESCO (Millero & Poisson 1981) equation of state — not linearised, not TEOS-10
- Frictional-geostrophic baroclinic velocity solve (algebraic, per-column, GPU-friendly)
- Barotropic streamfunction from a 2D elliptic vorticity equation with island constraints
- Vertical velocity diagnosed from ∇·**u** = 0 (rigid-lid approximation)
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

Oceananigans is used **only** as the tracer transport engine via
`PrescribedVelocityFields`. It does not solve momentum equations or maintain
a free surface — both are handled by Pelagos.

## Development phases

Development follows a strict sequential order validated against CLIMBER-X Fortran
reference output:

1. ✅ **EOS and pressure** — UNESCO EOS + hydrostatic integration
2. ✅ **Baroclinic velocity** — frictional-geostrophic algebraic solve
3. ✅ **Barotropic solver** — 2D elliptic ψ with island constraints
4. ✅ **Continuity** — diagnose w from ∇·u = 0
5. 🔲 **Oceananigans T/S model** — wire up PrescribedVelocityFields, closures
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
