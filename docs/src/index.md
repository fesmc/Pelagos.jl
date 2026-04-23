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

- UNESCO (Millero & Poisson 1981) equation of state — **not** linearised, **not** TEOS-10
- Frictional-geostrophic baroclinic velocity solve (algebraic, per-column)
- Barotropic streamfunction from an elliptic 2D vorticity equation with island constraints
- Vertical velocity diagnosed from ∇·**u** = 0 (rigid-lid approximation)
- Tracer (T, S) transport via [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl)
  with GM/Redi isopycnal diffusion and Bryan-Lewis diapycnal mixing
- Virtual salinity flux for freshwater forcing (rigid-lid consistent)
- Geothermal heat flux at the ocean floor

## Quick start

```julia
using Pelagos

# Compute seawater density at T=10°C, S=35 psu, p=50 bar
ρ = seawater_density(10.0, 35.0, 50.0)

# Build a simple Coriolis parameter array
f = coriolis_parameter.(collect(-80.0:5.0:80.0))
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
