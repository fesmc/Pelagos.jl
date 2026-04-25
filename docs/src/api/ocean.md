# Coupled Ocean Model

The `Ocean` module wraps the entire velocity-and-tracer pipeline. It is the
high-level entry point for time-stepping a Pelagos simulation.

`OceanModel` holds:

- the shared `ImmersedBoundaryGrid` (the single source of grid geometry),
- the Oceananigans `HydrostaticFreeSurfaceModel` that owns ``T`` and ``S`` and
  the prescribed velocity `Field`s,
- the pre-assembled sparse barotropic streamfunction solver, and
- the pressure `Field{Center, Center, Center}`.

`step_ocean!` advances the model one timestep using the pipeline:

1. [`compute_pressure!`](@ref Pelagos.Pressure.compute_pressure!) — hydrostatic
   ``p`` from ``T,S``.
2. [`solve_baroclinic!`](@ref Pelagos.Baroclinic.solve_baroclinic!) —
   frictional-geostrophic ``u, v`` on the C-grid.
3. Barotropic correction — solve ψ at T-points, average to corners, derive
   divergence-free face transports, replace the FG depth-mean with the
   ψ-derived barotropic flow.
4. [`diagnose_w!`](@ref Pelagos.Continuity.diagnose_w!) — vertical velocity
   from the rigid-lid continuity constraint.
5. `Oceananigans.time_step!` — advance ``T`` and ``S``.

The corner-ψ correction makes ``w_{\text{top}}`` discretely zero to
machine precision, so over one year of integration ``T`` and ``S`` drift
only at the diffusive rate (no spurious source/sink from a non-zero
rigid-lid residual).

```@autodocs
Modules = [Pelagos.Ocean]
```
