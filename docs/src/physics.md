# Physics Reference

This page summarises the governing equations implemented in Pelagos.jl. The
authoritative reference is **Appendix B of Willeit et al. (2022)**. Where the
paper and the CLIMBER-X Fortran source disagree, the Fortran is treated as
correct (see CLAUDE.md).

## Frictional-geostrophic momentum balance

At each interior grid column and vertical level, horizontal velocities satisfy:

```
-f·v + r_bc·u = -(1/ρ₀) ∂p/∂x + τˣ·δ(z=0)/h_top
 f·u + r_bc·v = -(1/ρ₀) ∂p/∂y + τʸ·δ(z=0)/h_top
```

Solving analytically for ``u`` and ``v``:

```math
u = \frac{r_\text{bc}\,F_x + f\,F_y}{r_\text{bc}^2 + f^2}, \qquad
v = \frac{r_\text{bc}\,F_y - f\,F_x}{r_\text{bc}^2 + f^2}
```

where ``F_x = -(1/\rho_0)\,\partial p/\partial x + \tau^x/h`` (top layer only)
and similarly for ``F_y``. The baroclinic friction coefficient is
``r_\text{bc} = 4\ \text{d}^{-1}``.

### Equatorial Coriolis floor

The frictional-geostrophic balance is singular at ``f = 0``. A floor is applied:

```math
f_\text{eff} = \text{sign}(f)\,\max(|f|,\, f_\text{min}), \quad f_\text{min} = 5\times10^{-6}\ \text{s}^{-1}
```

This floor is applied **identically** in both the baroclinic solve and the
barotropic vorticity equation.

## Barotropic streamfunction

The depth-integrated vorticity equation is:

```math
J\!\left(\psi,\, \frac{f}{H}\right) = \text{curl}\!\left(\frac{\boldsymbol{\tau}}{H}\right)
  - \frac{r_\text{bt}}{H}\,\nabla^2\psi + \text{baroclinic forcing}
```

``\psi`` is the barotropic streamfunction (m² s⁻¹), ``H`` is ocean depth,
``\boldsymbol{\tau}`` is wind stress, and ``r_\text{bt}`` is the spatially
variable barotropic friction. Near continental boundaries and in shallow water,
``r_\text{bt}`` is enhanced by a factor of 3.

Boundary conditions:
- ``\psi = 0`` on all connected continental boundaries.
- ``\psi = C_k`` (free constant) on each island ``k``.
- ``C_k`` is set by the island integral constraint: ``\oint (\partial\psi/\partial n)\,dl = 0``.

The system is assembled as a sparse linear system and solved with a direct solver.

## Continuity and vertical velocity

Vertical velocity is diagnosed by integrating upward from ``w = 0`` at the
ocean floor:

```math
w(k) = w(k+1) - \left(\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}\right)\Delta z_k
```

The rigid-lid constraint requires ``w_\text{surface} = 0`` to machine precision.

## Equation of state

The UNESCO equation of state of Millero & Poisson (1981) is implemented from
scratch. It expresses in-situ density as a rational function of temperature
(°C), salinity (psu), and pressure (bar):

```math
\rho(T, S, p) = \frac{\rho_0(T, S)}{1 - p/K(T, S, p)}
```

where ``\rho_0`` is the 1-atmosphere density and ``K`` is the secant bulk
modulus.

## Tracer diffusion

Lateral diffusion uses the Gent-McWilliams (GM) + Redi (isopycnal skew-symmetric)
scheme via Oceananigans' `IsopycnalSkewSymmetricDiffusivity`:

- Isopycnal diffusivity ``\kappa_\text{iso} = 1500\ \text{m}^2\ \text{s}^{-1}``
- GM coefficient ``\kappa_\text{GM} = 1500\ \text{m}^2\ \text{s}^{-1}``
- Gerdes et al. (1991) slope tapering at maximum slope ``10^{-3}``

Diapycnal diffusivity follows Bryan & Lewis (1979):

```math
\kappa_d(z) = \kappa_\text{bg} + (\kappa_\text{deep} - \kappa_\text{bg})\,\frac{2}{\pi}\,
   \arctan\!\left(\alpha\,(|z| - z_\text{ref})\right)
```

with ``\kappa_\text{bg} = 0.1\ \text{cm}^2\ \text{s}^{-1}`` near the surface
and ``\kappa_\text{deep} = 1.5\ \text{cm}^2\ \text{s}^{-1}`` at depth.

Convective adjustment uses `ConvectiveAdjustmentVerticalDiffusivity` with
``\kappa_\text{conv} = 100\ \text{m}^2\ \text{s}^{-1}``.

## Virtual salinity flux (rigid-lid freshwater forcing)

Under the rigid-lid approximation, freshwater forcing is implemented as a
virtual salinity flux:

```math
F_S = F_W \cdot S_\text{local} / h_\text{top}
```

where ``F_W`` is the net freshwater flux (m s⁻¹), ``S_\text{local}`` is the
surface salinity, and ``h_\text{top}`` is the top-layer thickness. A global
correction is applied each timestep to ensure the integral of ``F_S`` over the
ocean surface is zero.
