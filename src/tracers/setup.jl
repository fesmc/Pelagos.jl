# Oceananigans T/S model construction for Pelagos.jl.
#
# Architecture constraint: Oceananigans is used ONLY as a tracer transport
# engine.  Velocities (u, v, w) are prescribed each timestep by the Pelagos
# velocity solver; Oceananigans never solves momentum equations.
#
# The model is a HydrostaticFreeSurfaceModel with PrescribedVelocityFields.
# The ImplicitFreeSurface is a required API argument but is inactive when
# velocities are prescribed — Oceananigans skips the momentum step entirely.
#
# Buoyancy: SeawaterBuoyancy with a LinearEquationOfState is used solely so
# that the IsopycnalSkewSymmetricDiffusivity (GM/Redi) can compute isopycnal
# slopes.  The actual density used for the velocity solve comes from our own
# UNESCO EOS in src/eos/UNESCO.jl.
#
# Surface boundary conditions
# ───────────────────────────
# When `T_restore` / `S_restore` are provided, surface restoring is applied
# via FluxBoundaryCondition targeting the reference fields with timescale
# `tau_T` / `tau_S`.  The restoring flux is:
#
#   F_T = -dz_top · (T_sfc − T_ref) / τ_T        [°C m s⁻¹]
#   F_S = -dz_top · (S_sfc − S_ref) / τ_S        [psu m s⁻¹]
#
# where dz_top is the top layer thickness.  These have the sign convention of
# Oceananigans FluxBoundaryCondition (positive = into the ocean).

module TracerSetup

using Oceananigans
using Oceananigans.TurbulenceClosures
using Oceananigans.BoundaryConditions

using ..Diffusion: gm_redi_closure, diapycnal_closure
using ..Convection: convective_adjustment_closure
using ..Parameters: CP_OCN, RHO_0

export build_tracer_model, update_velocities!

# ── Model construction ─────────────────────────────────────────────────────────

"""
    build_tracer_model(grid;
        T_init    = nothing,
        S_init    = nothing,
        T_restore = nothing,
        S_restore = nothing,
        tau_T     = 15 * 86400.0,
        tau_S     = 15 * 86400.0,
    ) -> HydrostaticFreeSurfaceModel

Construct the Oceananigans tracer-only model on `grid`.

# Arguments
- `grid`      : Oceananigans grid (typically from `build_climberx_grid`)
- `T_init`    : initial temperature (nlon, nlat, nz) or nothing (→ 0)
- `S_init`    : initial salinity (nlon, nlat, nz) or nothing (→ 35)
- `T_restore` : 2-D surface temperature reference (nlon, nlat) for surface restoring
- `S_restore` : 2-D surface salinity reference (nlon, nlat) for surface restoring
- `tau_T`     : temperature restoring timescale, s  (default 15 days)
- `tau_S`     : salinity restoring timescale, s  (default 15 days)
"""
function build_tracer_model(grid;
                            T_init    = nothing,
                            S_init    = nothing,
                            T_restore = nothing,
                            S_restore = nothing,
                            tau_T     :: Float64 = 15.0 * 86400.0,
                            tau_S     :: Float64 = 15.0 * 86400.0)

    # ── Prescribed velocity fields (updated externally each timestep) ──────────
    u_field = Field{Face,   Center, Center}(grid)
    v_field = Field{Center, Face,   Center}(grid)
    w_field = Field{Center, Center, Face  }(grid)
    velocities = PrescribedVelocityFields(u=u_field, v=v_field, w=w_field)

    # ── Tracer closures ────────────────────────────────────────────────────────
    closure = (gm_redi_closure(), diapycnal_closure(), convective_adjustment_closure())

    # ── Buoyancy (linear EOS for GM/Redi slope computation only) ──────────────
    buoyancy = SeawaterBuoyancy(
        equation_of_state = LinearEquationOfState(
            thermal_expansion   = 2e-4,   # K⁻¹, typical for T≈10°C
            haline_contraction  = 8e-4,   # psu⁻¹
        )
    )

    # ── Surface boundary conditions ────────────────────────────────────────────
    T_bcs = _make_surface_bcs(T_restore, tau_T, grid, :T)
    S_bcs = _make_surface_bcs(S_restore, tau_S, grid, :S)

    boundary_conditions = (T = T_bcs, S = S_bcs)

    # ── Build model ────────────────────────────────────────────────────────────
    model = HydrostaticFreeSurfaceModel(;
        grid,
        velocities,
        closure,
        buoyancy,
        boundary_conditions,
        tracers            = (:T, :S),
        free_surface       = ImplicitFreeSurface(),
        coriolis           = nothing,
        momentum_advection = nothing,
        tracer_advection   = WENO(),
    )

    # ── Initial conditions ─────────────────────────────────────────────────────
    if !isnothing(T_init)
        set!(model.tracers.T, _flip_to_oceananigans(T_init))
    end
    if !isnothing(S_init)
        set!(model.tracers.S, _flip_to_oceananigans(S_init))
    end

    return model
end

# ── Velocity update ────────────────────────────────────────────────────────────

"""
    update_velocities!(model, u_C, v_C, w)

Copy C-grid velocity arrays into the Oceananigans PrescribedVelocityFields.

# Arguments
- `u_C` : (Nλ+1, Nφ,   Nz)   zonal velocity on east faces, m s⁻¹
- `v_C` : (Nλ,   Nφ+1, Nz)   meridional velocity on north faces, m s⁻¹
- `w`   : (Nλ,   Nφ,   Nz+1) vertical velocity on z-faces, m s⁻¹
"""
function update_velocities!(model,
                            u_C::AbstractArray{Float64,3},
                            v_C::AbstractArray{Float64,3},
                            w  ::AbstractArray{Float64,3})
    interior(model.velocities.u) .= u_C
    interior(model.velocities.v) .= v_C
    interior(model.velocities.w) .= w
    return nothing
end

# ── Internal helpers ───────────────────────────────────────────────────────────

# CLIMBER-X arrays are (nlon, nlat, nz) with k=1 at the deepest layer, which
# matches Oceananigans' convention — no reordering is needed.
_flip_to_oceananigans(arr) = arr

# Build a FieldBoundaryConditions for surface restoring, or return defaults.
function _make_surface_bcs(restore_ref, tau, grid, tracer_sym)
    isnothing(restore_ref) && return FieldBoundaryConditions()

    # dz_top = thickness of the top (surface) layer
    Nz  = grid.Nz
    dz_top = @inbounds grid.underlying_grid.Δzᵃᵃᶜ[Nz]

    params = (ref=restore_ref, tau=tau, dz=dz_top)

    @inline function restoring_flux(i, j, grid, clock, fields, p)
        tracer_field = tracer_sym == :T ? fields.T : fields.S
        sfc_val = @inbounds tracer_field[i, j, grid.Nz]
        ref_val = @inbounds p.ref[i, j]
        return -p.dz * (sfc_val - ref_val) / p.tau
    end

    top_bc = FluxBoundaryCondition(restoring_flux;
        discrete_form = true,
        parameters    = params,
    )
    return FieldBoundaryConditions(top = top_bc)
end

end # module TracerSetup
