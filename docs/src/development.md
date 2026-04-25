# Development Guide

## Running the tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Or from the Julia REPL:

```julia
using Pkg
Pkg.test("Pelagos")
```

Tests are organised by development phase:

```
test/
├── eos/             # Phase 1: UNESCO EOS and pressure
├── velocity/        # Phases 2–4: velocity solver
│   ├── test_islands.jl
│   ├── test_baroclinic.jl
│   ├── test_barotropic.jl
│   └── test_continuity.jl
├── tracers/         # Phase 5: diffusion closures, forcing
├── integration/     # Phase 6: spinup (requires fixtures)
└── fixtures/        # Reference NetCDF output (not in git)
```

## Fixture-based validation

Integration tests require NetCDF reference output from CLIMBER-X. See
`test/fixtures/README.md` for how to generate and place these files.

## Adding a new module

1. Create the source file in the appropriate `src/` subdirectory.
2. Add `include` and `using` statements in `src/Pelagos.jl`.
3. Export public functions in `src/Pelagos.jl`.
4. Create a corresponding test file in `test/`.
5. Include the test file in `test/runtests.jl`.
6. Add documentation in `docs/src/api/`.

## Coding conventions

- All 3D arrays: `(nlon, nlat, nz)` — longitude varies fastest.
- No magic numbers: all constants in `src/parameters.jl`.
- Type-stable functions: verify with `@code_warntype`.
- No module-level mutable globals.
- Default precision: `Float64` throughout.
- No explicit horizontal loops in the baroclinic solve — use broadcasting.
- See `CLAUDE.md` for the full list of conventions.

## Referencing the Fortran source

When implementing a physics module, read the corresponding CLIMBER-X Fortran
file first:

| Julia module | Fortran reference |
|---|---|
| `src/eos/UNESCO.jl` | — (pure Millero & Poisson 1981) |
| `src/eos/pressure.jl` | `ocn_pressure.f90` |
| `src/velocity/baroclinic.jl` | `ocn_baroclinic.f90` |
| `src/velocity/barotropic.jl` | `ocn_barotropic.f90` |
| `src/velocity/islands.jl` | `ocn_islands.f90` |
| `src/forcing/surface.jl` | `ocn_forcing.f90` |
| `src/grid/bathymetry.jl` | `ocn_grid.f90` |

Fortran files are in the CLIMBER-X repository at `src/ocn/`.

## Validation tolerances

| Variable | Tolerance |
|---|---|
| Density ρ | < 1×10⁻⁴ kg m⁻³ |
| Baroclinic u, v | < 1×10⁻⁴ m s⁻¹ RMS |
| Barotropic ψ | < 0.5 Sv |
| ``w_{\text{top}}`` (rigid-lid residual) | < 1×10⁻¹² m s⁻¹ (machine precision) |
| ``w`` interior | < 1×10⁻⁷ m s⁻¹ RMS |
| T (year 100) | < 0.5 °C RMS |
| S (year 100) | < 0.1 psu RMS |
| AMOC at 26°N (year 100) | within 3 Sv |

The rigid-lid residual is now constrained at machine precision because the
corner-ψ barotropic correction makes the depth-integrated transport
discretely divergence-free. See [Architecture](@ref) for details.
