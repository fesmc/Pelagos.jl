# Oceananigans T/S model construction for Pelagos.jl.
#
# Wires the custom velocity solver output (u, v, w as plain arrays) into
# Oceananigans via PrescribedVelocityFields, then constructs a HydrostaticModel
# (or NonhydrostaticModel with forcings disabled) that handles only tracer
# advection-diffusion.
#
# Architecture constraint: Oceananigans MUST NOT solve momentum equations.
# We use PrescribedVelocityFields so Oceananigans sees prescribed u,v,w and
# evolves only T and S.

module TracerSetup

using Oceananigans
using Oceananigans.Models: HydrostaticFreeSurfaceModel
using Oceananigans.TurbulenceClosures
using Oceananigans.Grids: RectilinearGrid, LatitudeLongitudeGrid

using ..Diffusion: gm_redi_closure, diapycnal_closure
using ..Convection: convective_adjustment_closure

export build_tracer_model, update_velocities!

"""
    build_tracer_model(grid; T_init=nothing, S_init=nothing) -> model

Construct the Oceananigans tracer-only model on the given `grid`.

The model uses PrescribedVelocityFields (initialised to zero) for u,v,w.
These are updated externally each timestep via `update_velocities!`.

Closures: GM/Redi + diapycnal Bryan-Lewis + convective adjustment.

# Arguments
- `grid`   : an Oceananigans grid (LatitudeLongitudeGrid or RectilinearGrid)
- `T_init` : initial temperature field (nothing → zero)
- `S_init` : initial salinity field (nothing → zero)
"""
function build_tracer_model(grid;
                            T_init = nothing,
                            S_init = nothing)
    # Prescribed (externally updated) velocity fields
    u_field = Field{Face,  Center, Center}(grid)
    v_field = Field{Center, Face,  Center}(grid)
    w_field = Field{Center, Center, Face}(grid)

    velocities = PrescribedVelocityFields(u = u_field,
                                          v = v_field,
                                          w = w_field)

    # Combined closures: GM/Redi + diapycnal + convection
    closure = (gm_redi_closure(), diapycnal_closure(), convective_adjustment_closure())

    model = HydrostaticFreeSurfaceModel(;
        grid,
        velocities,
        closure,
        tracers   = (:T, :S),
        buoyancy  = nothing,        # no buoyancy forcing through Oceananigans
        free_surface = nothing,     # rigid lid: no Oceananigans free surface
        momentum_advection = nothing,
    )

    # Set initial conditions if provided
    T_init !== nothing && set!(model.tracers.T, T_init)
    S_init !== nothing && set!(model.tracers.S, S_init)

    return model
end

"""
    update_velocities!(model, u_arr, v_arr, w_arr)

Copy velocity arrays (nlon, nlat, nz) from the custom solver into the
Oceananigans PrescribedVelocityFields backing the tracer model.
"""
function update_velocities!(model,
                            u_arr::AbstractArray{Float64,3},
                            v_arr::AbstractArray{Float64,3},
                            w_arr::AbstractArray{Float64,3})
    interior(model.velocities.u) .= u_arr
    interior(model.velocities.v) .= v_arr
    # w has nz+1 faces; map accordingly
    interior(model.velocities.w)[:,:,1:end-1] .= w_arr[:,:,1:end-1]
    return nothing
end

end # module TracerSetup
