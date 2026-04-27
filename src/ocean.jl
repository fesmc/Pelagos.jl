# Full ocean model — wires the Phase 1-4 velocity solver with the Oceananigans
# T/S tracer transport engine.
#
# All state lives in Oceananigans Fields on a shared ImmersedBoundaryGrid.
# The same grid is used by every module (pressure, baroclinic, barotropic,
# continuity, tracer model), so there is no k-convention flip and no separate
# GridParams struct.
#
# Oceananigans k-convention: k=1 = deepest layer, k=Nz = shallowest layer.
# Velocity locations: u on (Face, Center, Center),
#                    v on (Center, Face, Center),
#                    w on (Center, Center, Face).
# Pressure p on (Center, Center, Center) — T-points.
#
# Pipeline each timestep:
#   1. compute_pressure!   — p from T, S
#   2. solve_baroclinic!   — u, v from p, wind stress
#   3. Barotropic correction: solve ψ, strip FG depth-mean, add ψ-derived flow
#   4. diagnose_w!         — w from ∇·u = 0
#   5. update velocities in tracer model, time_step! T, S

module Ocean

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: AbstractGrid, Center, Face
using Oceananigans.Fields: Field, interior
using Oceananigans.Operators: Δxᶜᶜᶜ, Δxᶜᶠᶜ, Δyᶜᶜᶜ, Δzᵃᵃᶜ
using Oceananigans.Grids: inactive_node, peripheral_node

using ..Pressure: compute_pressure!
using ..Baroclinic: solve_baroclinic!, coriolis_parameter
using ..Barotropic: BarotropicSolver, build_barotropic_solver, solve_barotropic!,
                     compute_jebar_forcing!
using ..Continuity: diagnose_w!
using ..TracerSetup: build_tracer_model, update_velocities!
using ..GridSetup: build_climberx_grid
using ..Restarts: load_climberx_restart

export OceanModel, build_ocean_model, step_ocean!

# ── Model struct ───────────────────────────────────────────────────────────────

"""
    OceanModel

Holds all state for one timestep of the coupled velocity + tracer ocean model.

Fields:
- `grid`         : shared ImmersedBoundaryGrid (the single source of grid geometry)
- `tracer_model` : Oceananigans HydrostaticFreeSurfaceModel (owns T, S and their
                   prescribed velocity fields u, v, w)
- `bar_solver`   : pre-assembled sparse barotropic streamfunction solver
- `p`            : hydrostatic pressure field (Center, Center, Center), bar
"""
mutable struct OceanModel
    grid         :: AbstractGrid
    tracer_model :: Any                            # HydrostaticFreeSurfaceModel
    bar_solver   :: BarotropicSolver
    p            :: Field{Center, Center, Center}
    ubar_jbar    :: Matrix{Float64}                # JEBAR forcing scratch (Nx, Ny)
end

# ── Constructor ────────────────────────────────────────────────────────────────

"""
    build_ocean_model(restart_dir; T_restore, S_restore, tau_T, tau_S)

Build an OceanModel initialised from a CLIMBER-X restart directory.
The grid is constructed from the restart file; initial T and S are loaded
into the Oceananigans tracer fields via `set!`.
"""
function build_ocean_model(restart_dir ::String;
                           T_restore              = nothing,
                           S_restore              = nothing,
                           tau_T      ::Float64   = 15.0 * 86400.0,
                           tau_S      ::Float64   = 15.0 * 86400.0)

    restart_file = joinpath(restart_dir, "ocn_restart.nc")

    # Build the Oceananigans grid from the restart file (zw, k1_pot, f_ocn)
    grid = build_climberx_grid(CPU(); restart_file)

    # Load initial T, S (Oceananigans k-convention: k=1=deepest already)
    r = load_climberx_restart(restart_dir)

    # Build tracer model (HydrostaticFreeSurfaceModel with PrescribedVelocityFields)
    tracer_model = build_tracer_model(grid;
        T_init    = r.T,
        S_init    = r.S,
        T_restore,
        S_restore,
        tau_T,
        tau_S,
    )

    # Barotropic solver assembled once from grid geometry
    bar_solver = build_barotropic_solver(grid)

    # Pressure field (T-point, same grid as tracers)
    p = Field{Center, Center, Center}(grid)

    # JEBAR forcing scratch buffer (recomputed each step from current T, S)
    ubar_jbar = zeros(Float64, grid.Nx, grid.Ny)

    return OceanModel(grid, tracer_model, bar_solver, p, ubar_jbar)
end

# ── Single timestep ────────────────────────────────────────────────────────────

"""
    step_ocean!(model, tau_x, tau_y, dt)

Advance the ocean model by one timestep `dt` seconds.

`tau_x`, `tau_y` are (Nx, Ny) matrices of wind stress in N m⁻² on T-points.
"""
function step_ocean!(m     ::OceanModel,
                     tau_x ::Matrix{Float64},
                     tau_y ::Matrix{Float64},
                     dt    ::Float64)
    grid = m.grid
    Nx   = grid.Nx
    Ny   = grid.Ny
    Nz   = grid.Nz

    # Convenience handles into the tracer model's velocity fields
    u = m.tracer_model.velocities.u   # Field{Face,   Center, Center}
    v = m.tracer_model.velocities.v   # Field{Center, Face,   Center}
    w = m.tracer_model.velocities.w   # Field{Center, Center, Face  }
    T = m.tracer_model.tracers.T      # Field{Center, Center, Center}
    S = m.tracer_model.tracers.S      # Field{Center, Center, Center}

    # ── 1. Hydrostatic pressure ────────────────────────────────────────────────
    # compute_pressure! fills halos on p, including copying pressure into
    # immersed cells so horizontal gradients at the seafloor stay small.
    compute_pressure!(m.p, T, S, grid)

    # ── 2. Baroclinic frictional-geostrophic velocities ───────────────────────
    solve_baroclinic!(u, v, m.p, tau_x, tau_y, grid)

    # ── 3. Barotropic correction ───────────────────────────────────────────────
    # Goal: replace the FG depth-mean (which has non-zero depth-integrated
    # divergence on the C-grid) with a ψ-derived barotropic flow that is
    # EXACTLY divergence-free in the discrete sense used by diagnose_w!.
    #
    # Key discrete fact: if ψ lives at cell CORNERS (F-points) and the C-grid
    # face transports are
    #     U(i-½, j) = −(ψ_c[i, j+1] − ψ_c[i, j])       (flux, m³/s)
    #     V(i, j-½) =  (ψ_c[i+1, j] − ψ_c[i,   j])
    # then the finite-volume divergence (δx Ax·u + δy Ay·v) at T-cell (i,j)
    # telescopes to zero identically — independent of H, metric factors, or
    # how ψ_c was obtained.  Proof: expand the 4 flux differences and cancel.
    #
    # The existing barotropic solver produces ψ at T-points (it solves the
    # vorticity equation on the T-grid for physics reasons).  We average T-point
    # ψ to corners with a 4-point bilinear average; ψ_corner = 0 on the south
    # and north domain walls.  This loses a half-cell of accuracy in the ψ
    # field but guarantees exact discrete divergence freedom of the barotropic
    # transport.  A future refactor can solve the vorticity equation directly
    # on the F-grid.

    ocean_mask, H, f_arr, dx, dy = _extract_barotropic_geometry(grid)

    # JEBAR forcing — implementation complete (compute_jebar_forcing! uses
    # face-bottom-summed integration with a level-mean density anomaly,
    # giving |J(1/H, E)| ≈ wind-curl magnitude on PI ρ, T, S — see
    # docstring).  Wiring it in below currently destabilises the no-wind
    # tracer integration after ~120 days: the ~30 Sv barotropic flow
    # JEBAR adds, combined with the existing FG baroclinic velocities,
    # overshoots the first-order upwind tracer advection.  Re-enable
    # together with a higher-order monotone tracer advection or a
    # smaller dt.
    #
    # compute_jebar_forcing!(m.ubar_jbar, T, S, grid)
    # psi_T = solve_barotropic!(m.bar_solver, tau_x, tau_y, H, f_arr, dx, dy, ocean_mask;
    #                           ubar_jbar = m.ubar_jbar)
    psi_T = solve_barotropic!(m.bar_solver, tau_x, tau_y, H, f_arr, dx, dy, ocean_mask)

    # Map ψ_T to corners.  ψ_c[i, j] is the corner at position (i−½, j−½),
    # with i=1..Nx (periodic in x) and j=1..Ny+1.  Boundary corners (j=1 and
    # j=Ny+1) are zero — the south and north domain walls.
    #
    # At coastal corners (any of the 4 surrounding T-cells is land), ψ_c is set
    # to 0 — the continent's streamfunction value.  Without this, averaging
    # mixes psi_T=0 (land) with non-zero ocean values, breaking the boundary
    # condition "ψ = constant on continent" and giving non-zero discrete
    # divergence at coastal T-cells.
    psi_c = zeros(Nx, Ny+1)
    @inbounds for j in 2:Ny, i in 1:Nx
        im1 = mod1(i-1, Nx)
        # Are all 4 surrounding T-cells ocean?
        all_ocean = ocean_mask[im1, j-1] && ocean_mask[i, j-1] &&
                    ocean_mask[im1, j  ] && ocean_mask[i, j  ]
        if all_ocean
            psi_c[i, j] = 0.25 * (psi_T[im1, j-1] + psi_T[i, j-1] +
                                  psi_T[im1, j  ] + psi_T[i, j  ])
        end
        # Coastal/boundary corner: psi_c stays 0 (continent constant)
    end

    # Column depths at U-faces and V-faces (sum of Δz over active layers).
    H_U = zeros(Nx, Ny)
    H_V = zeros(Nx, Ny+1)
    u_d = interior(u)
    v_d = interior(v)
    u_FG_mean = zeros(Nx, Ny)       # depth-mean of FG u at each U-face
    v_FG_mean = zeros(Nx, Ny+1)     # depth-mean of FG v at each V-face

    @inbounds for k in 1:Nz
        for j in 1:Ny, i in 1:Nx
            # U-face (i-½, j) is active if not immersed at this level
            if !peripheral_node(i, j, k, grid, Face(), Center(), Center())
                dzk = Δzᵃᵃᶜ(i, j, k, grid)
                H_U[i, j]       += dzk
                u_FG_mean[i, j] += u_d[i, j, k] * dzk
            end
        end
        for j in 2:Ny, i in 1:Nx
            # V-face (i, j-½) is active if not immersed
            if !peripheral_node(i, j, k, grid, Center(), Face(), Center())
                dzk = Δzᵃᵃᶜ(i, j, k, grid)
                H_V[i, j]       += dzk
                v_FG_mean[i, j] += v_d[i, j, k] * dzk
            end
        end
    end
    @inbounds for j in 1:Ny,   i in 1:Nx; H_U[i,j] > 0 && (u_FG_mean[i,j] /= H_U[i,j]); end
    @inbounds for j in 2:Ny,   i in 1:Nx; H_V[i,j] > 0 && (v_FG_mean[i,j] /= H_V[i,j]); end

    # ψ-derived barotropic velocities at C-grid faces.
    # At U-face (i−½, j): u_bt = −(ψ_c[i, j+1] − ψ_c[i, j]) / (Δy · H_U)
    # At V-face (i, j−½): v_bt = +(ψ_c[i+1, j] − ψ_c[i, j]) / (Δx_Vface · H_V)
    u_bt = zeros(Nx, Ny)
    v_bt = zeros(Nx, Ny+1)
    @inbounds for j in 1:Ny, i in 1:Nx
        H_U[i,j] > 0.0 || continue
        u_bt[i, j] = -(psi_c[i, j+1] - psi_c[i, j]) / (dy * H_U[i, j])
    end
    @inbounds for j in 2:Ny, i in 1:Nx
        H_V[i,j] > 0.0 || continue
        ip1 = mod1(i+1, Nx)
        # Δx at V-face j uses the V-face latitude (carries cosφ metric)
        Δx_Vface = Δxᶜᶠᶜ(i, j, 1, grid)
        v_bt[i, j] = (psi_c[ip1, j] - psi_c[i, j]) / (Δx_Vface * H_V[i, j])
    end

    # Apply correction at every level: strip FG depth-mean, add ψ-derived.
    # u_total[k] = u_FG[k] − u_FG_mean + u_bt
    # v_total[k] = v_FG[k] − v_FG_mean + v_bt
    @inbounds for k in 1:Nz
        for j in 1:Ny, i in 1:Nx
            peripheral_node(i, j, k, grid, Face(), Center(), Center()) && continue
            u_d[i, j, k] += (u_bt[i, j] - u_FG_mean[i, j])
        end
        for j in 2:Ny, i in 1:Nx
            peripheral_node(i, j, k, grid, Center(), Face(), Center()) && continue
            v_d[i, j, k] += (v_bt[i, j] - v_FG_mean[i, j])
        end
    end

    fill_halo_regions!(u)
    fill_halo_regions!(v)

    # ── 4. Vertical velocity from continuity ───────────────────────────────────
    diagnose_w!(w, u, v, grid)

    # ── 5. Advance T, S ────────────────────────────────────────────────────────
    time_step!(m.tracer_model, dt)

    # Clamp salinity after advection to prevent WENO undershoots from
    # accumulating and causing UNESCO sqrt(S) DomainError.
    interior(S) .= max.(interior(S), 0.0)

    return nothing
end

# ── Internal helpers ───────────────────────────────────────────────────────────

function _extract_barotropic_geometry(grid::AbstractGrid)
    Nx = grid.Nx;  Ny = grid.Ny
    bh = grid.immersed_boundary.bottom_height   # (Nx, Ny), m, negative = ocean

    ocean_mask = Matrix{Bool}(undef, Nx, Ny)
    H          = zeros(Float64, Nx, Ny)
    dx         = zeros(Float64, Ny)

    dy = Δyᶜᶜᶜ(1, 1, 1, grid)

    @inbounds for j in 1:Ny
        dx[j] = Δxᶜᶜᶜ(1, j, 1, grid)
        for i in 1:Nx
            is_ocean       = bh[i,j] < 0.0
            ocean_mask[i,j] = is_ocean
            H[i,j]          = is_ocean ? -bh[i,j] : 0.0
        end
    end

    f_arr = [begin
                 φ = grid.φᵃᶜᵃ[j]
                 fj = 2.0 * (2π/86164.0) * sind(φ)
                 s  = fj >= 0.0 ? 1.0 : -1.0
                 s * max(abs(fj), 5e-6)
             end for j in 1:Ny]

    return ocean_mask, H, f_arr, dx, dy
end

end # module Ocean
