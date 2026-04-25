# Architecture

## Overview

```
┌──────────────────────────────────────────────────────────────────────┐
│                              Pelagos.jl                              │
│                                                                      │
│  Single ImmersedBoundaryGrid + Oceananigans Fields shared by all:    │
│                                                                      │
│  ┌──────────────────────────┐         ┌─────────────────────────┐    │
│  │   Custom velocity solver │ u,v,w   │  Oceananigans tracer    │    │
│  │   (operates on Fields)   │ ──────► │  model (T, S advection) │    │
│  │                          │         │                         │    │
│  │   1. UNESCO ρ(T,S,p)     │         │  PrescribedVelocityFlds │    │
│  │   2. ∫ρg dz → p (Field)  │         │  GM/Redi (Iso skew)     │    │
│  │   3. FG baroclinic u,v   │         │  Bryan-Lewis κ_d(z)     │    │
│  │      ∂xᶠᶜᶜ p, ∂yᶜᶠᶜ p    │         │  Convective adjustment  │    │
│  │   4. Barotropic ψ solve  │         │  UpwindBiased(order=1)  │    │
│  │      + corner-ψ correct  │         │                         │    │
│  │   5. w = ∇·u = 0 (flux)  │         │                         │    │
│  └──────────────────────────┘         └─────────────────────────┘    │
└──────────────────────────────────────────────────────────────────────┘
```

The full pipeline is wrapped by `OceanModel` and `step_ocean!` in
`src/ocean.jl`.

## Fundamental constraint

**GOLDSTEIN is not a Navier-Stokes model.** Oceananigans' built-in models
solve prognostic momentum equations; GOLDSTEIN diagnoses velocities
algebraically. The solution is to use Oceananigans' tracer engine via
`PrescribedVelocityFields`, supplying externally computed ``u``, ``v``, ``w``
every timestep — but to do so on the **same `Field`s, same grid, same halos,
same immersed-boundary masking** that Oceananigans uses internally.

## Field-and-grid architecture

All physics state is stored in Oceananigans `Field`s on a shared
`ImmersedBoundaryGrid` (built once in `GridSetup.build_climberx_grid`):

| Field | Location | Owner |
|---|---|---|
| ``T``, ``S`` | `(Center, Center, Center)` | tracer model |
| ``u`` | `(Face, Center, Center)` | tracer model (prescribed) |
| ``v`` | `(Center, Face, Center)` | tracer model (prescribed) |
| ``w`` | `(Center, Center, Face)` | tracer model (prescribed) |
| ``p`` | `(Center, Center, Center)` | `OceanModel` |

Consequences:

- **No `(nlon, nlat, nz)` arrays as a parallel state.** Everything is a
  `Field`; finite differences use `Oceananigans.Operators` (`∂xᶠᶜᶜ`,
  `∂yᶜᶠᶜ`, `Δzᵃᵃᶜ`, `Axᶠᶜᶜ`, …) that carry the correct spherical metrics.
- **No k-convention flip.** All modules use Oceananigans' k=1=deepest
  convention. Restart data from CLIMBER-X is loaded directly without
  reordering.
- **Halo regions handle boundaries.** `fill_halo_regions!(p)` takes care
  of periodic-x and Bounded-y conditions; the FG solver and `diagnose_w!`
  index `Field`s rather than raw `interior(...)` arrays so halo cells
  (e.g. across the dateline) are accessible.
- **Immersed-boundary masking.** `inactive_node` and `peripheral_node`
  from `Oceananigans.Grids` decide which T-cells / faces are active.

## Velocity-face masking: `peripheral_node`, not `inactive_node`

For the FG balance to give meaningful velocities at coastal columns with
**stepped bathymetry**, a U-face must be closed (``u = 0``) when *either*
neighbouring T-cell is below bathymetry — the `peripheral_node` predicate.
Using only `inactive_node` (true when *both* neighbours are below
bathymetry) leaves U-faces "open" alongside a step, where the pressure
gradient mixes pressures from incompatible depths and produces unphysical
flow that contaminates the depth-integrated divergence.

## Corner-ψ barotropic correction

The barotropic streamfunction ``\psi`` is solved on T-points by the
existing physics (sparse vorticity equation with island constraints).
For the **discrete** depth-integrated continuity to be exactly zero on a
C-grid, ``\psi`` is averaged from T-points to F-points (cell corners),
with ``\psi_c = 0`` enforced at any corner adjacent to a continent. The
C-grid face transports

```math
U_{i-\tfrac12,\,j} = -(\psi_c[i,j+1] - \psi_c[i,j]),\qquad
V_{i,\,j-\tfrac12} = +(\psi_c[i+1,j] - \psi_c[i,j])
```

then make the four flux differences at every T-cell telescope to zero
algebraically — independent of ``H``, metric factors, or how ``\psi_c``
was obtained. The barotropic correction strips the FG depth-mean from
each face and adds these ``\psi_c``-derived velocities; the resulting
3-D field has machine-precision rigid-lid residual ``w_{\text{top}}``.

(A future refactor can solve the vorticity equation directly on the
F-grid; the current corner-averaging loses a half-cell of accuracy in
``\psi`` itself but keeps the discrete divergence exact.)

## Continuity: custom flux-form, not `div_xyᶜᶜᶜ`

`Oceananigans.Operators.div_xyᶜᶜᶜ` on an `ImmersedBoundaryGrid` applies
a `conditional_δ` that **zeros the entire divergence** at any T-cell
adjacent to a closed face. For our rigid-lid integration we want the
opposite behaviour: open-face fluxes still contribute, while closed-face
fluxes vanish naturally because ``u = v = 0`` there. `Continuity.diagnose_w!`
therefore uses a custom flux-form divergence on the underlying grid
(`Axᶠᶜᶜ`, `Ayᶜᶠᶜ`) that bypasses the conditional masking.

## Oceananigans integration

Oceananigans owns:
- ``T`` and ``S`` advection-diffusion
- GM/Redi + diapycnal + convective closures
- The grid (`LatitudeLongitudeGrid` + `GridFittedBottom`)
- All `Field`s, halo management, and finite-difference operators

Oceananigans does **not** own:
- ``u``, ``v``, ``w`` (prescribed each step from the FG solve + ψ correction)
- Pressure (diagnosed from ``\rho(T,S,p)`` via UNESCO)
- The free surface (rigid-lid model — no ``\eta``)

Tracer advection uses `UpwindBiased(order=1)` to match GOLDSTEIN's
upwind scheme and to avoid the WENO5 undershoot that would drive
``S < 0`` and trigger the UNESCO ``\sqrt{S}`` domain error.

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
| `Ocean` | `src/ocean.jl` | Phase 5 |

## Development phases

Development follows a strict sequential order; each phase must pass its tests
before the next begins.

1. ✅ **EOS and pressure** — UNESCO EOS + hydrostatic integration
2. ✅ **Baroclinic velocity** — frictional-geostrophic algebraic solve
3. ✅ **Barotropic solver** — 2D elliptic ψ with island constraints
4. ✅ **Continuity** — diagnose ``w`` from ``\nabla\cdot\mathbf{u} = 0``
5. ✅ **Oceananigans T/S model** — coupled `OceanModel` runs stably for
   ≥1 year; ``w_{\text{top}}`` at machine precision
6. 🔲 **Integration validation** — 100-year spinup vs. CLIMBER-X fixtures

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
