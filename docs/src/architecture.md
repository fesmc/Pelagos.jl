# Architecture

## Overview

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

## Fundamental constraint

**GOLDSTEIN is not a Navier-Stokes model.** Oceananigans solves prognostic
momentum equations; GOLDSTEIN diagnoses velocities algebraically. These are
architecturally incompatible. The solution is to use Oceananigans **only** as
a tracer transport engine via `PrescribedVelocityFields`, supplying externally
computed ``u``, ``v``, ``w`` every timestep.

## Oceananigans integration

Oceananigans owns:
- ``T`` and ``S`` advection-diffusion
- GM/Redi + diapycnal + convective closures
- Grid and Field definitions

Oceananigans does **not** own:
- ``u``, ``v``, ``w`` (prescribed externally)
- Pressure (diagnosed from ``\rho(T,S,p)``)
- The free surface (rigid-lid model)

## Module map

| Module | File | Phase |
|---|---|---|
| `Parameters` | `src/parameters.jl` | All |
| `UNESCO` | `src/eos/UNESCO.jl` | Phase 1 |
| `Pressure` | `src/eos/pressure.jl` | Phase 1 |
| `Islands` | `src/velocity/islands.jl` | Phase 3 |
| `Baroclinic` | `src/velocity/baroclinic.jl` | Phase 2 |
| `Barotropic` | `src/velocity/barotropic.jl` | Phase 3 |
| `Continuity` | `src/velocity/continuity.jl` | Phase 4 |
| `Bathymetry` | `src/grid/bathymetry.jl` | Phase 1 |
| `Masks` | `src/grid/masks.jl` | Phase 1 |
| `Surface` | `src/forcing/surface.jl` | Phase 5 |
| `Geothermal` | `src/forcing/geothermal.jl` | Phase 5 |
| `Diffusion` | `src/tracers/diffusion.jl` | Phase 5 |
| `Convection` | `src/tracers/convection.jl` | Phase 5 |
| `TracerSetup` | `src/tracers/setup.jl` | Phase 5 |
| `Output` | `src/io/output.jl` | Phase 5 |
| `Restarts` | `src/io/restarts.jl` | Phase 5 |

## Development phases

Development follows a strict sequential order; each phase must pass its tests
before the next begins.

1. **EOS and pressure** — UNESCO EOS + hydrostatic integration
2. **Baroclinic velocity** — frictional-geostrophic algebraic solve
3. **Barotropic solver** — 2D elliptic ψ with island constraints
4. **Continuity** — diagnose ``w`` from ``\nabla\cdot\mathbf{u} = 0``
5. **Oceananigans T/S model** — wire up `PrescribedVelocityFields`, closures
6. **Integration validation** — 100-year spinup vs. CLIMBER-X fixtures

## Known physics traps

See the full list in `CLAUDE.md`. Key items:

- **Rigid lid ≠ free surface**: the rigid-lid pressure diagnostic is different
  from Oceananigans' sea surface height.
- **Virtual salinity flux**: applying real freshwater flux causes long-term salinity drift.
- **Equatorial f-floor consistency**: ``f_\text{min}`` must be identical in
  the baroclinic solve and the barotropic vorticity equation.
- **Spatially variable barotropic friction**: near boundaries, ``r_\text{bt}``
  is multiplied by 3.
- **Tracer/velocity mask distinction**: Bering and Davis Straits are closed for
  barotropic flow but open for tracer exchange.
