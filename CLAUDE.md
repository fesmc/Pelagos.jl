# CLAUDE.md — Pelagos.jl

This file is the authoritative guide for all development work on this repository.
Read it in full before writing any code. It takes precedence over your general
knowledge of Oceananigans, Julia ocean modelling conventions, or anything inferred
from partial context.

---

## What this project is

A Julia reimplementation of the GOLDSTEIN frictional-geostrophic ocean model,
targeting functional equivalence with the GOLDSTEIN component of CLIMBER-X v1.0
(Willeit et al., GMD 2022, https://doi.org/10.5194/gmd-15-5905-2022).
The associated pdf is in the repository here:
willeit-et-al_2022_the-earth-system-model-climber-x-v1-0-part-1.pdf

The reference Fortran source is the CLIMBER-X repository. The authoritative
equation reference is **Appendix B of the CLIMBER-X paper**. Where the paper and
the Fortran disagree, the Fortran is correct.

The goal is a Julia package `Pelagos.jl` that:
1. Reproduces CLIMBER-X GOLDSTEIN ocean state (T, S, overturning streamfunction,
   barotropic streamfunction) to within acceptable tolerance when driven by the
   same surface forcing and bathymetry.
2. Uses Oceananigans.jl as its tracer transport engine only, as well as for 
   grid and Fields definitions, and other helper functions.
3. Is written in idiomatic Julia with clear module boundaries and unit tests.

---

## Architecture — read this before touching any code

### The fundamental constraint

**GOLDSTEIN is not a Navier-Stokes model. Do not use Oceananigans' momentum solver.**

GOLDSTEIN diagnoses horizontal velocities algebraically from a frictional-geostrophic
balance at every timestep. There is no prognostic momentum equation. Oceananigans
solves prognostic momentum equations. These are incompatible at the core.

The chosen architecture is:

```
┌─────────────────────────────────────────────────────────┐
│                     Pelagos.jl                          │
│                                                         │
│  ┌──────────────┐    u,v,w     ┌─────────────────────┐  │
│  │  VelocitySolver │ ────────► │  Oceananigans        │  │
│  │  (custom Julia) │           │  PrescribedVelocity  │  │
│  │                 │           │  Fields + T,S model  │  │
│  │  1. EOS + ∂p/∂z │           │                      │  │
│  │  2. Baroclinic  │           │  Owns: T, S          │  │
│  │     FG solve    │           │  advection-diffusion │  │
│  │  3. Barotropic  │           │  GM/Redi, convection │  │
│  │     ψ solve     │           │                      │  │
│  │  4. w from ∇·u=0│           │  Does NOT own:       │  │
│  └──────────────┘            │  u, v, w, p          │  │
│                               └─────────────────────┘  │
└─────────────────────────────────────────────────────────┘
```

Oceananigans is used **only** via `PrescribedVelocityFields`. It owns T and S.
It does not own u, v, w, or pressure. The velocity solver is entirely custom Julia
and lives in `src/velocity/`.

### What Oceananigans must not be asked to do

- Do not use `HydrostaticFreeSurfaceModel` for the velocity field.
- Do not use `NonhydrostaticModel`.
- Do not use Oceananigans' Coriolis closure on momentum.
- Do not use Oceananigans' pressure solver for the ocean pressure.
- Do not use Oceananigans' built-in `SplitExplicitFreeSurface` or
  `ImplicitFreeSurface`.

If you find yourself reaching for any of these, stop and reconsider.

### Rigid lid and virtual salinity flux

GOLDSTEIN uses a **rigid-lid** approximation. There is no free surface.
Freshwater forcing is implemented as a **virtual salinity flux**:

```
FS = FW * S_local / H_toplayer
```

where FW is the net freshwater flux (m/s), S_local is the local surface salinity,
and H_toplayer is the top layer thickness. This must be applied as the salinity
surface boundary condition. Do **not** apply a real freshwater volume flux. This
is a known difference from free-surface models that must be preserved for
conservation.

A global correction is applied each timestep to ensure the globally integrated
virtual salinity flux is zero (see Appendix B, Section B4).

---

## Module structure (to be determined if it is the best way)

```
src/
├── Pelagos.jl              # Package root, re-exports public API
├── eos/
│   ├── UNESCO.jl           # Millero & Poisson (1981) equation of state
│   └── pressure.jl         # Hydrostatic pressure integration
├── velocity/
│   ├── baroclinic.jl       # Frictional-geostrophic baroclinic solve
│   ├── barotropic.jl       # 2D barotropic streamfunction solver
│   ├── continuity.jl       # Diagnose w from ∇·u = 0
│   └── islands.jl          # Island integral constraints
├── tracers/
│   ├── setup.jl            # Oceananigans T/S model construction
│   ├── diffusion.jl        # GM/Redi + diapycnal diffusivity profiles
│   └── convection.jl       # Convective adjustment
├── forcing/
│   ├── surface.jl          # Wind stress, heat flux, virtual salinity flux
│   └── geothermal.jl       # Bottom boundary condition (Lucazeau 2019)
├── grid/
│   ├── bathymetry.jl       # Bathymetry loading and smoothing
│   └── masks.jl            # Land mask, island detection
└── io/
    ├── output.jl           # NetCDF output
    └── restarts.jl         # Checkpoint read/write

test/
├── eos/                    # EOS unit tests vs. tabulated values
├── velocity/               # Velocity solver tests vs. reference snapshots
│   ├── baroclinic/
│   ├── barotropic/
│   └── continuity/
├── tracers/                # Tracer transport tests
├── integration/            # Multi-century spinup vs. reference NetCDF
└── fixtures/               # Reference NetCDF output from CLIMBER-X Fortran
    ├── README.md           # Documents how fixtures were generated
    ├── snapshot_year0010.nc
    ├── snapshot_year0100.nc
    └── snapshot_year1000.nc
```

Each module has a corresponding test directory. **Do not move on to the next
module until the tests for the current one pass against the fixtures.**

---

## Physics — key equations and decisions

### Frictional-geostrophic baroclinic momentum balance

At each interior grid column and level:

```
-f·v + r_bc·u = -(1/ρ₀) ∂p/∂x + τˣ·δ(z=0)/h_top
 f·u + r_bc·v = -(1/ρ₀) ∂p/∂y + τʸ·δ(z=0)/h_top
```

where:
- `f` is the Coriolis parameter (with floor, see below)
- `r_bc = 4 d⁻¹` is the globally uniform baroclinic friction coefficient
- `p` is the hydrostatic pressure computed from ρ(T,S,p) via the EOS
- Wind stress `τ` is applied only at the top layer, divided by layer thickness

Solving the 2×2 system analytically:

```
u = (r_bc·Fx + f·Fy) / (r_bc² + f²)
v = (r_bc·Fy - f·Fx) / (r_bc² + f²)
```

where `Fx = -(1/ρ₀)∂p/∂x + τˣ/h` and similarly for `Fy`.

This solve is **local** (per column, per level) and embarrassingly parallel.
Vectorise over the horizontal grid; do not loop in Julia.

### Equatorial Coriolis floor

The frictional-geostrophic balance is singular at f = 0. Apply a floor:

```julia
f_eff = sign(f) * max(abs(f), f_min)
```

where `f_min = 5e-6 s⁻¹`. This must be applied **consistently** in both the
baroclinic solve (above) and the barotropic vorticity equation. Do not apply it
anywhere else (e.g. not in the GM scheme's f-dependence if that is ever added).

### Barotropic streamfunction

The depth-integrated flow is governed by:

```
J(ψ, f/H) = curl(τ/H) - r_bt·∇²ψ/H + (baroclinic_forcing)
```

This is an elliptic equation for the barotropic streamfunction ψ. It must be
solved on the full 2D ocean domain (land cells excluded) with:

- `ψ = 0` on all connected land boundaries
- One free constant per island (non-simply-connected domain regions)
- Island constants determined by the integral constraint:
  `∮ (∂ψ/∂n) dl = 0` around each island perimeter

The barotropic friction `r_bt` is spatially variable: it is multiplied by a
factor of 3 near continental boundaries and shallow bathymetric features (see
Appendix B of the paper and the Fortran source `src/ocean/ocn_barotropic.f90`).

**Islands with non-zero barotropic flow** (controlled by namelist):
- Drake Passage: enabled by default
- Indonesian Throughflow: enabled by default
- Bering Strait: **closed** for barotropic flow (but open for tracers)
- Davis Strait: **closed** for barotropic flow (but open for tracers)

Island detection must be performed automatically from the bathymetry mask. See
`src/grid/masks.jl`. Do not hardcode island indices.

The elliptic solve should use a sparse direct solver (`SparseArrays` +
`LinearSolve.jl` with a suitable backend). At 5°×5° the system is ~2500 unknowns —
direct solvers are appropriate; do not use iterative solvers unless benchmarking
shows a clear need.

The barotropic streamfunction must be recomputed **every timestep**, not cached.

### Vertical velocity from continuity

After u and v are known at all levels, w is diagnosed by integrating the continuity
equation upward from the bottom (w=0 at bathymetry):

```
w(k) = w(k-1) - (∂u/∂x + ∂v/∂y)·Δz(k)
```

with metric factors for the spherical grid. The top-level w should be zero to
machine precision (rigid lid); if it is not, there is a bug in the barotropic
solve. Add an assertion in debug mode.

### Tracer diffusion

Use Oceananigans' `IsopycnalSkewSymmetricDiffusivity` for the combined
GM + Redi scheme. Parameters:

- Isopycnal diffusivity κ_iso = 1500 m² s⁻¹
- GM coefficient κ_gm = 1500 m² s⁻¹ (equal to κ_iso by default)
- Isoneutral slope tapering following Gerdes et al. (1991):
  maximum slope = 1×10⁻³

Diapycnal diffusivity follows Bryan & Lewis (1979):

```
κ_d(z) = κ_bg + (κ_deep - κ_bg) * (2/π) * arctan(α·(|z| - z_ref))
```

with κ_bg = 0.1 cm² s⁻¹ at surface, κ_deep = 1.5 cm² s⁻¹ at ocean floor.
Implement as a `VerticalScalarDiffusivity` with depth-dependent profile.

### Convective adjustment

Use Oceananigans' `ConvectiveAdjustmentVerticalDiffusivity` with a large
convective diffusivity (e.g. κ_conv = 100 m² s⁻¹). The Rahmstorf (1993) scheme
used in the original GOLDSTEIN is effectively equivalent at the coarse resolution
targeted here.

### Equation of state

Implement the UNESCO equation of state of Millero & Poisson (1981). Do **not**
use Oceananigans' built-in `LinearEquationOfState` or `TEOS10`. The UNESCO EOS
must be implemented from scratch in `src/eos/UNESCO.jl`. Validate against the
check values in Millero & Poisson (1981) Table 5 before use.

### Bathymetry smoothing

The bathymetry must be smoothed by convolving with a 4-point nearest-neighbour
kernel before use (see Section 2.2 of the paper). This is applied once during
grid initialisation, not at each timestep. Additionally, a limitation is applied
to the topographic slope entering the Coriolis term in the barotropic equation;
see the Fortran source `src/ocean/ocn_grid.f90` for the exact implementation.

---

## Reference Fortran files

The ocean-relevant CLIMBER-X Fortran source files are located in this subdirectory of 
the CLIMBER-X model repository:

https://github.com/cxesmc/climber-x/src/ocn/

When implementing a Julia module, **read the corresponding Fortran file first**.
Where the Fortran and the paper differ, the Fortran is authoritative. Note
discrepancies in comments in the Julia source.

---

## Validation and testing

### Philosophy

Every module must be independently testable against known-correct output before
integration. Do not proceed to the next phase if tests are failing.

### Reference fixtures

Pre-computed NetCDF reference output from the original CLIMBER-X Fortran code
lives in `test/fixtures/`. See `test/fixtures/README.md` for details on how
they were generated and what variables are included. Key variables:

- `ocn_temp`, `ocn_salt` — 3D T and S fields
- `ocn_u`, `ocn_v`, `ocn_w` — 3D velocity fields
- `ocn_psi_bt` — barotropic streamfunction (2D)
- `ocn_psi_moc` — meridional overturning streamfunction (zonal mean)
- `ocn_pressure` — 3D hydrostatic pressure
- `ocn_rho` — 3D density

### Test tolerances

| Variable | Tolerance | Notes |
|---|---|---|
| Density ρ | < 1×10⁻⁴ kg m⁻³ | vs. Millero & Poisson table values |
| Baroclinic u, v | < 1×10⁻⁴ m s⁻¹ RMS | vs. snapshot fixture |
| Barotropic ψ | < 0.5 Sv | vs. snapshot fixture |
| w | < 1×10⁻⁷ m s⁻¹ RMS | vs. snapshot fixture |
| T (year 100) | < 0.5 °C RMS | vs. spinup fixture |
| S (year 100) | < 0.1 psu RMS | vs. spinup fixture |
| AMOC at 26°N (year 100) | within 3 Sv | vs. spinup fixture |

Tighter is better. These are minimum acceptance thresholds, not targets.

### Required tests per module

**`src/eos/UNESCO.jl`**
- Reproduce Table 5 check values from Millero & Poisson (1981)
- Reproduce `ocn_rho` from `snapshot_year0010.nc` within tolerance

**`src/velocity/baroclinic.jl`**
- Given `ocn_pressure` and surface wind stress from fixture, reproduce `ocn_u`, `ocn_v`

**`src/velocity/barotropic.jl`**
- Given depth-integrated forcing from fixture, reproduce `ocn_psi_bt`
- Verify that island integrals are satisfied (net barotropic transport through
  Drake Passage matches fixture; Bering Strait transport is zero)

**`src/velocity/continuity.jl`**
- Given `ocn_u`, `ocn_v` from fixture, reproduce `ocn_w`
- Verify top-level w < 1×10⁻¹² m s⁻¹ (rigid-lid constraint)

**Integration test**
- 10-year spinup from rest with present-day forcing
- 100-year spinup from present-day initial conditions
- Compare to `snapshot_year0010.nc` and `snapshot_year0100.nc`

---

## Julia conventions for this project

### Column-major and vectorisation

All 3D arrays are stored `(nlon, nlat, nz)` — longitude varies fastest in memory
(column-major, consistent with Oceananigans' default). The baroclinic solve must
be vectorised over `(nlon, nlat)` at each level using broadcasting; do not write
explicit loops over horizontal indices in Julia.

Use Oceananigans.Fields and related operators as much as possible. Allow for different architectures (CPU/GPU).

### No magic numbers

All physical constants and tuning parameters must be defined in
`src/parameters.jl` with units in comments. Do not embed literal values in
physics functions.

```julia
const R_BC     = 4.0 / 86400.0   # s⁻¹, baroclinic friction coefficient
const R_BT_FAC = 3.0              # amplification near boundaries
const F_MIN    = 5e-6             # s⁻¹, Equatorial Coriolis floor
const K_ISO    = 1500.0           # m² s⁻¹, isopycnal diffusivity
const K_GM     = 1500.0           # m² s⁻¹, GM coefficient
const SLOPE_MAX = 1e-3            # Gerdes et al. slope taper threshold
```

### Type stability

All physics functions must be type-stable (verify with `@code_warntype`). Use
concrete-typed structs for parameter collections. Avoid `Any`-typed fields.

### No global mutable state

Do not use module-level mutable globals for model state. All state lives in
structs passed explicitly.

### Float precision

Default to `Float64` throughout. If single precision is ever explored, it must
be via a type parameter, not by changing constants.

---

## Dependencies

Pinned in `Project.toml`. Do not add dependencies without discussion.

| Package | Purpose |
|---|---|
| `Oceananigans.jl` | Tracer transport engine (T, S only) |
| `LinearSolve.jl` | Barotropic elliptic solver |
| `SparseArrays` | Standard library, sparse matrix assembly |
| `NCDatasets.jl` | NetCDF I/O for forcing, output, fixtures |
| `Interpolations.jl` | Depth-profile interpolation (Bryan-Lewis κ) |
| `Test` | Standard library, unit tests |

---

## Known physics traps — read before implementing

**1. The rigid-lid / free-surface distinction is not cosmetic.**
The rigid-lid pressure that GOLDSTEIN implies (used for sea ice in CLIMBER-X)
is diagnosed by integrating density above 1500 m depth, not by solving a surface
pressure equation. If you ever need sea surface height, use this diagnostic, not
a free-surface height from Oceananigans.

**2. Virtual salinity flux is not a minor detail.**
Applying real freshwater flux instead of virtual salinity flux will cause salinity
drift that compounds over millennium-scale integrations. The global correction
term (ensuring net annual flux is zero) is also mandatory for conservation. See
Appendix B4 and `src/ocean/ocn_forcing.f90`.

**3. The barotropic friction enhancement is spatially variable.**
Near continental boundaries and shallow features, r_bt is multiplied by 3. The
definition of "near boundary" follows the Fortran implementation in
`ocn_barotropic.f90`. Do not apply uniform barotropic friction.

**4. Island integrals change when the land-sea mask changes.**
If the bathymetry or land-sea mask is ever updated at runtime (e.g. for
paleo-simulations with evolving land-sea distribution), the island topology
must be redetected and the barotropic solver must be re-assembled. The sparse
matrix for the elliptic solve is not constant across mask changes.

**5. Pressure gradient discretisation on the sphere.**
The pressure gradient in the frictional-geostrophic solve must use the same
metric factors and grid staggering as the Fortran. In particular, ∂p/∂x on the
sphere is (1/a·cosφ)·∂p/∂λ. Check the Fortran `ocn_pressure.f90` carefully for
the staggering convention (pressure at cell centres, velocities at cell faces).

**6. The Equatorial f floor must be consistent.**
The same `f_min = 5e-6 s⁻¹` floor must be used in both the baroclinic algebraic
solve and in the barotropic vorticity equation. Using it in one but not the other
will produce a barotropic/baroclinic inconsistency that may not fail loudly but
will cause incorrect overturning.

**7. Tracer exchange through closed straits.**
Bering Strait and Davis Strait are closed for barotropic flow (ψ = const across
them) but open for tracer exchange (T and S can diffuse through). This is
controlled by separate masks for velocity and tracers. Do not use the same mask
for both.

---

## Development phases

Work must proceed in this order. Do not start a phase until the previous phase's
tests pass.

**Phase 1 — EOS and pressure** (`src/eos/`)
Implement UNESCO EOS and hydrostatic pressure integration. No Oceananigans
dependency. Tests: Millero & Poisson table values; reproduce `ocn_rho` and
`ocn_pressure` from fixture.

**Phase 2 — Baroclinic velocity** (`src/velocity/baroclinic.jl`)
Given pressure and wind stress, produce u and v. No Oceananigans dependency.
Tests: reproduce `ocn_u`, `ocn_v` from fixture.

**Phase 3 — Barotropic solver** (`src/velocity/barotropic.jl`, `islands.jl`)
2D elliptic solve for ψ with island constraints. No Oceananigans dependency.
Tests: reproduce `ocn_psi_bt` from fixture; verify Drake Passage and Indonesian
Throughflow transports; verify Bering and Davis closure.

**Phase 4 — Continuity** (`src/velocity/continuity.jl`)
Diagnose w from ∇·u = 0. Tests: reproduce `ocn_w`; verify rigid-lid constraint
at surface.

**Phase 5 — Oceananigans T/S model** (`src/tracers/`)
Wire velocity solver output into `PrescribedVelocityFields`. Implement GM/Redi,
diapycnal diffusion, convection, surface/bottom boundary conditions. Tests:
10-year spinup within tolerance.

**Phase 6 — Integration validation**
100-year spinup. Compare T, S, AMOC to fixtures within tolerances defined above.

---

## Questions that require human input before proceeding

If you encounter any of the following, stop and ask rather than making a
decision:

- The Fortran source for a given module is ambiguous or appears inconsistent with
  the paper equations.
- A test is failing and the cause is unclear after reading both the Fortran and
  the paper.
- A new Julia dependency seems necessary.
- The island topology of a paleo bathymetry differs significantly from present
  day in a way that breaks the current island detection logic.
- Any change to the architecture described above seems necessary or beneficial.
- Tolerances in the validation tests are not met and it is not obvious why.
