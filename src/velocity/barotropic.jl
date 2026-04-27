# Barotropic streamfunction solver for GOLDSTEIN ocean model.
#
# Solves the elliptic vorticity equation for the barotropic streamfunction ψ:
#   J(ψ, f/H) = curl(τ/H) - r_bt·∇²ψ/H + baroclinic_forcing
#
# where J is the Jacobian operator, H is ocean depth, r_bt is the (spatially
# variable) barotropic friction coefficient, τ is the wind stress vector.
#
# Boundary conditions:
#   ψ = 0      on all connected-continent boundaries
#   ψ = C_k    (free constant) on each island k
#   C_k determined by: ∮ (∂ψ/∂n) dl = 0  around island k perimeter
#
# The system is assembled as a sparse linear system and solved with LinearSolve.jl.
# At 5°×5° the system is ~2500 unknowns; direct solvers are appropriate.
#
# Barotropic friction is enhanced by R_BT_FAC near continental boundaries and
# in shallow water (see CLIMBER-X ocn_barotropic.f90 for exact criterion).
#
# Reference: Appendix B of Willeit et al. (2022) and ocn_barotropic.f90.

module Barotropic

using SparseArrays
using Oceananigans
using Oceananigans.Grids: AbstractGrid
using Oceananigans.Operators: Δxᶜᶜᶜ, Δyᶜᶜᶜ, Δzᵃᵃᶜ
using Oceananigans.Grids: inactive_node
using Oceananigans.Fields: Field, Center, interior

using ..Parameters: R_BT_BASE, R_BT_FAC, F_MIN, RHO_0, G, R_EARTH
using ..UNESCO: seawater_density
using ..Islands: IslandInfo, detect_islands

export BarotropicSolver, build_barotropic_solver, compute_rbt, solve_barotropic!
export compute_jebar_forcing, compute_jebar_forcing!

"""
    BarotropicSolver

Pre-assembled sparse system for the barotropic streamfunction.  Assembled once
per land-sea mask and reused every timestep.

Fields:
- `A`         : sparse coefficient matrix
- `rhs`       : right-hand side vector (updated each timestep)
- `psi_vec`   : solution vector (nlon*nlat unknowns, land cells zeroed)
- `cell_index`: (nlon, nlat) → linear index into the ocean-cell DOF ordering
- `islands`   : detected island list
"""
mutable struct BarotropicSolver
    A         :: SparseMatrixCSC{Float64, Int}
    rhs       :: Vector{Float64}
    psi_vec   :: Vector{Float64}
    cell_index:: Matrix{Int}          # 0 = land; >0 = DOF index
    islands   :: Vector{IslandInfo}
    n_dof     :: Int
end

"""
    build_barotropic_solver(ocean_mask, f, H, dx, dy, lat_deg) -> BarotropicSolver

Assemble the sparse linear system for the barotropic streamfunction.

# Arguments
- `ocean_mask`: Bool (nlon, nlat)
- `f`         : Coriolis (nlat,), s⁻¹, with F_MIN floor applied
- `H`         : ocean depth (nlon, nlat), m, positive downward
- `dx`        : zonal grid spacing (nlat,), m
- `dy`        : meridional grid spacing, m (uniform assumed)
- `r_bt`      : barotropic friction field (nlon, nlat), s⁻¹
"""
function build_barotropic_solver(ocean_mask::Matrix{Bool},
                                 f         ::Vector{Float64},
                                 H         ::Matrix{Float64},
                                 dx        ::Vector{Float64},
                                 dy        ::Float64,
                                 r_bt      ::Matrix{Float64})::BarotropicSolver
    nlon, nlat = size(ocean_mask)
    islands    = detect_islands(ocean_mask)

    # Assign DOF indices to ocean cells + island DOFs
    cell_index = zeros(Int, nlon, nlat)
    n_ocean    = 0
    for j in 1:nlat, i in 1:nlon
        ocean_mask[i, j] || continue
        n_ocean += 1
        cell_index[i, j] = n_ocean
    end
    n_islands = length(islands)
    n_dof     = n_ocean + n_islands

    # Sparse triplets (I, J, V)
    Is = Int[]; Js = Int[]; Vs = Float64[]
    rhs = zeros(Float64, n_dof)

    # Finite-volume stencil at each ocean T-cell.
    #
    # The vorticity equation
    #     J(ψ, f/H) − ∇·((r/H) ∇ψ) = (1/ρ₀) curl(τ/H)
    # is discretised in flux form: each T-face carries its own (r/H)·∇ψ
    # flux and its own (f/H) value used in the Jacobian central difference.
    # H_face = min(H_left, H_right) — the CLIMBER-X convention so a face
    # cannot be deeper than the shallower of its bounding cells.  Land-
    # adjacent faces are closed (zero flux, zero coefficient).
    #
    # f varies only in latitude, so:
    #   ∂(f/H)/∂x at T-cell centre  =  f(j) · (1/H_E − 1/H_W) / (2 Δx)
    #   ∂(f/H)/∂y at T-cell centre  =  (f_N/H_N − f_S/H_S) / (2 Δy)
    # with f_N = ½(f(j) + f(jp1)), f_S = ½(f(j) + f(jm1)) at the v-faces.
    for j in 1:nlat, i in 1:nlon
        ocean_mask[i, j] || continue
        row = cell_index[i, j]
        Hi  = max(H[i, j], 1.0)
        ri  = r_bt[i, j]
        dxj = dx[j]

        ip1 = mod1(i + 1, nlon); im1 = mod1(i - 1, nlon)
        jp1 = min(j + 1, nlat);  jm1 = max(j - 1, 1)

        # ── Face geometry ──────────────────────────────────────────────────
        # GOLDSTEIN/CLIMBER-X convention: H_face = min(H_left, H_right) at
        # open faces.  At closed faces (neighbour is land, ψ_land = 0) we
        # project the cell-centre value onto the face — this preserves the
        # diagonal contribution that enforces the Dirichlet BC ψ = 0 on
        # continents, while giving zero off-diagonal coefficient.
        open_E = ocean_mask[ip1, j];  open_W = ocean_mask[im1, j]
        open_N = ocean_mask[i, jp1] && j < nlat
        open_S = ocean_mask[i, jm1] && j > 1

        H_E = open_E ? min(Hi, H[ip1, j]) : Hi
        H_W = open_W ? min(Hi, H[im1, j]) : Hi
        H_N = open_N ? min(Hi, H[i, jp1]) : Hi
        H_S = open_S ? min(Hi, H[i, jm1]) : Hi

        r_E = open_E ? 0.5 * (ri + r_bt[ip1, j]) : ri
        r_W = open_W ? 0.5 * (ri + r_bt[im1, j]) : ri
        r_N = open_N ? 0.5 * (ri + r_bt[i, jp1]) : ri
        r_S = open_S ? 0.5 * (ri + r_bt[i, jm1]) : ri

        # ── Friction (flux form) ───────────────────────────────────────────
        # ∇·((r/H) ∇ψ) at T-cell, with each face contributing (r/H)_face/Δ²
        Cfric_E = r_E / (H_E * dxj^2)
        Cfric_W = r_W / (H_W * dxj^2)
        Cfric_N = r_N / (H_N * dy^2)
        Cfric_S = r_S / (H_S * dy^2)
        diag    = -(Cfric_E + Cfric_W + Cfric_N + Cfric_S)

        # ── Jacobian (face-staggered f/H) ──────────────────────────────────
        # ∂(f/H)/∂x at T-centre  =  f(j) · (1/H_E − 1/H_W)/(2 Δx)
        # ∂(f/H)/∂y at T-centre  =  (f_N/H_N − f_S/H_S)/(2 Δy)
        invH_E = H_E > 0 ? 1.0/H_E : 0.0
        invH_W = H_W > 0 ? 1.0/H_W : 0.0
        invH_N = H_N > 0 ? 1.0/H_N : 0.0
        invH_S = H_S > 0 ? 1.0/H_S : 0.0
        f_N = 0.5 * (f[j] + f[jp1])
        f_S = 0.5 * (f[j] + f[jm1])
        dfH_dx = f[j] * (invH_E - invH_W) / (2.0 * dxj)
        dfH_dy = (f_N * invH_N - f_S * invH_S) / (2.0 * dy)

        # J(ψ, f/H) ≈ (∂ψ/∂y)·dfH_dx − (∂ψ/∂x)·dfH_dy
        # Central differences of ψ → coefficients on N/S/E/W neighbours.
        coeff_J_NS = dfH_dx / (2.0 * dy)     # multiplies (ψ_N − ψ_S)
        coeff_J_EW = dfH_dy / (2.0 * dxj)    # multiplies −(ψ_E − ψ_W)

        # ── Diagonal ───────────────────────────────────────────────────────
        # All 4 faces closed = isolated ocean cell: pin ψ = 0.  (Avoids
        # singular row when the system is solved.)
        if !(open_E || open_W || open_N || open_S)
            push!(Is, row); push!(Js, row); push!(Vs, 1.0)
            continue
        end
        push!(Is, row); push!(Js, row); push!(Vs, diag)

        # ── Off-diagonals: friction (+) and Jacobian (±) combined ─────────
        if open_E
            _add_neighbour!(Is, Js, Vs, row, cell_index, ip1, j, ocean_mask,
                             Cfric_E - coeff_J_EW)
        end
        if open_W
            _add_neighbour!(Is, Js, Vs, row, cell_index, im1, j, ocean_mask,
                             Cfric_W + coeff_J_EW)
        end
        if open_N
            _add_neighbour!(Is, Js, Vs, row, cell_index, i, jp1, ocean_mask,
                             Cfric_N + coeff_J_NS)
        end
        if open_S
            _add_neighbour!(Is, Js, Vs, row, cell_index, i, jm1, ocean_mask,
                             Cfric_S - coeff_J_NS)
        end
    end

    # Island constraint rows: ∮ (∂ψ/∂n) dl = 0
    # These are additional equations (one per island) with ψ_island as the unknown.
    for (k, isl) in enumerate(islands)
        island_dof = n_ocean + k
        # Perimeter integral:  sum over perimeter cells of (ψ_perim - C_k) / dist = 0
        # → sum(ψ_perim) / N - C_k = 0
        row = island_dof
        n_perim = length(isl.perimeter)
        if n_perim > 0
            coeff_each = 1.0 / n_perim
            for (ip, jp) in isl.perimeter
                col = cell_index[ip, jp]
                col > 0 || continue
                push!(Is, row); push!(Js, col); push!(Vs, coeff_each)
            end
        end
        push!(Is, row); push!(Js, island_dof); push!(Vs, -1.0)
        rhs[row] = 0.0
    end

    A = sparse(Is, Js, Vs, n_dof, n_dof)
    psi_vec = zeros(Float64, n_dof)
    return BarotropicSolver(A, rhs, psi_vec, cell_index, islands, n_dof)
end

@inline function _add_neighbour!(Is, Js, Vs, row, cell_index, ni, nj, ocean_mask, coeff)
    ocean_mask[ni, nj] || return
    idx = cell_index[ni, nj]
    idx > 0 || return
    push!(Is, row); push!(Js, idx); push!(Vs, coeff)
end

"""
    compute_rbt(ocean_mask, H, dx, dy) -> Matrix{Float64}

Compute spatially variable barotropic friction: R_BT_BASE × R_BT_FAC near
continental boundaries and in shallow water (H < 500 m), R_BT_BASE elsewhere.
"""
function compute_rbt(ocean_mask::Matrix{Bool},
                     H         ::Matrix{Float64},
                     dx        ::Vector{Float64},
                     dy        ::Float64)::Matrix{Float64}
    nlon, nlat = size(ocean_mask)
    r_bt = fill(R_BT_BASE, nlon, nlat)

    for j in 1:nlat, i in 1:nlon
        ocean_mask[i, j] || continue
        # Shallow-water criterion
        if H[i, j] < 500.0
            r_bt[i, j] = R_BT_BASE * R_BT_FAC
            continue
        end
        # Continental-boundary criterion: any 4-connected land neighbour
        for (di, dj) in ((-1,0),(1,0),(0,-1),(0,1))
            ni = mod1(i + di, nlon); nj = j + dj
            (1 ≤ nj ≤ nlat) || continue
            if !ocean_mask[ni, nj]
                r_bt[i, j] = R_BT_BASE * R_BT_FAC
                break
            end
        end
    end
    return r_bt
end

"""
    solve_barotropic!(solver, tau_x, tau_y, H, f, dx, dy, ocean_mask) -> Matrix{Float64}

Solve for the barotropic streamfunction ψ in m³ s⁻¹ (1 Sv = 10⁶ m³ s⁻¹).
Returns ψ on the full (nlon, nlat) grid (land cells = 0).

Wind stress `tau_x`, `tau_y` in N m⁻²; the assembler divides by `RHO_0` so
the equation is dimensionally consistent.
"""
function solve_barotropic!(solver     ::BarotropicSolver,
                           tau_x      ::Matrix{Float64},
                           tau_y      ::Matrix{Float64},
                           H          ::Matrix{Float64},
                           f          ::Vector{Float64},
                           dx         ::Vector{Float64},
                           dy         ::Float64,
                           ocean_mask ::Matrix{Bool};
                           ubar_jbar  ::Union{Nothing, Matrix{Float64}} = nothing
                          )::Matrix{Float64}
    nlon, nlat = size(ocean_mask)
    rhs   = solver.rhs
    ci    = solver.cell_index

    fill!(rhs, 0.0)

    # Assemble RHS: (1/ρ₀) · [curl(τ/H) + JEBAR].
    # The barotropic vorticity equation has units s⁻²:
    #   J(ψ, f/H)  −  ∇·((r_bt/H) ∇ψ)  =  (1/ρ₀) · curl(τ/H) − (1/ρ₀) · J(E, 1/H)
    # where E = ∫ p_bc dz is the depth-integrated baroclinic pressure
    # (Pa·m).  Both forcings have units kg·m⁻³·s⁻² (= Pa·m⁻²); /ρ₀ brings
    # them to s⁻² to match the LHS.
    #
    # JEBAR sign: the Jacobian J(1/H, E) = -J(E, 1/H), so we ADD ubar_jbar
    # to the wind-stress term — `ubar_jbar` itself is already J(1/H, E).
    for j in 1:nlat, i in 1:nlon
        ocean_mask[i, j] || continue
        row = ci[i, j]
        Hi  = max(H[i, j], 1.0)
        dxj = dx[j]

        ip1 = mod1(i+1,nlon); im1 = mod1(i-1,nlon)
        jp1 = min(j+1,nlat);  jm1 = max(j-1,1)

        # curl(τ/H) = ∂(τy/H)/∂x − ∂(τx/H)/∂y
        dtauY_dx = (tau_y[ip1,j]/max(H[ip1,j],1.0) - tau_y[im1,j]/max(H[im1,j],1.0)) / (2.0*dxj)
        dtauX_dy = (tau_x[i,jp1]/max(H[i,jp1],1.0) - tau_x[i,jm1]/max(H[i,jm1],1.0)) / (2.0*dy)
        wind_term = dtauY_dx - dtauX_dy

        jbar_term = isnothing(ubar_jbar) ? 0.0 : ubar_jbar[i, j]

        rhs[row] = (wind_term + jbar_term) / RHO_0
    end

    # Solve A·ψ = rhs using Julia's built-in sparse direct solver
    solver.A[end, end] != 0.0 || (solver.A[1,1] += 1e-10)  # ensure non-singular
    psi_vec = solver.A \ rhs

    # Unpack solution onto 2D grid
    psi = zeros(Float64, nlon, nlat)
    for j in 1:nlat, i in 1:nlon
        ocean_mask[i, j] || continue
        idx = ci[i, j]
        idx > 0 && (psi[i, j] = psi_vec[idx])
    end
    return psi
end

# ── JEBAR forcing ────────────────────────────────────────────────────────────

"""
    compute_jebar_forcing(T, S, grid) -> Matrix{Float64}

Compute the JEBAR (Joint Effect of Baroclinicity And Relief) forcing for the
barotropic-vorticity RHS, in units of kg m⁻³ s⁻² (= Pa m⁻², equivalent to
the wind-stress curl term *before* division by ρ₀).

The CLIMBER-X formulation (jbar.f90):

  bp(i, j, k) = − ∫_{z_bottom}^{z_k} ρ g dz        (Pa, zero at the bottom)

For each face F between two T-cells, integrate `bp · dz` from a *common
face-bottom level* `k_face = max(k_bot_left, k_bot_right)` — using the
deeper of the two cells' active ranges so any sub-shelf contributions
that would differ between the two cells are excluded.  Then form

  E_F     = Σ_{k≥k_face} ½(bp_left(k) + bp_right(k)) dz(k)        (Pa·m)
  invH_F  = 1 / Σ_{k≥k_face} dz(k)                                 (m⁻¹)

The T-cell JEBAR is the central difference on those face values:

  ubar_jbar[i,j] = ∂(1/H)/∂x · ∂E/∂y − ∂(1/H)/∂y · ∂E/∂x

with ∂(1/H)/∂x = (invH_E − invH_W)/Δx (and similarly for y), and likewise
for ∂E.  At ocean–land faces the face-bottom level lies above any active
layer, so both E_F and invH_F vanish identically — the closed-boundary
treatment falls out of the staggering rather than a special case.

Periodic in λ; bounded north/south.
"""
function compute_jebar_forcing(T    ::Field{Center, Center, Center},
                               S    ::Field{Center, Center, Center},
                               grid ::AbstractGrid)::Matrix{Float64}
    out = zeros(Float64, grid.Nx, grid.Ny)
    compute_jebar_forcing!(out, T, S, grid)
    return out
end

function compute_jebar_forcing!(ubar_jbar::Matrix{Float64},
                                T    ::Field{Center, Center, Center},
                                S    ::Field{Center, Center, Center},
                                grid ::AbstractGrid)
    Nx = grid.Nx; Ny = grid.Ny; Nz = grid.Nz
    Td = interior(T); Sd = interior(S)

    # ── Pre-step: density anomaly ρ' = ρ − ρ̄(k) ─────────────────────────
    # Subtract the horizontal-mean density at each level.  J(1/H, ρ̄·H²) is
    # mathematically zero (ρ̄ depends only on z), but its discrete residual
    # on a coarse grid swamps the physical baroclinic signal by ~4 orders
    # of magnitude.  Removing the level-mean removes that noise floor.
    # k_bot[i,j] = deepest active level (Nz+1 sentinel for fully-land).
    ρ_anom = zeros(Float64, Nx, Ny, Nz)
    k_bot  = fill(Nz + 1, Nx, Ny)
    @inbounds for j in 1:Ny, i in 1:Nx
        for k in 1:Nz
            inactive_node(i, j, k, grid, Center(), Center(), Center()) && continue
            k_bot[i, j] = k; break
        end
    end
    @inbounds for k in 1:Nz
        sum_ρ = 0.0; n = 0
        for j in 1:Ny, i in 1:Nx
            inactive_node(i, j, k, grid, Center(), Center(), Center()) && continue
            ρij = seawater_density(Td[i,j,k], max(Sd[i,j,k], 0.0), 0.0)
            sum_ρ += ρij; n += 1
        end
        ρ_mean_k = n > 0 ? sum_ρ / n : 0.0
        for j in 1:Ny, i in 1:Nx
            inactive_node(i, j, k, grid, Center(), Center(), Center()) && continue
            ρij = seawater_density(Td[i,j,k], max(Sd[i,j,k], 0.0), 0.0)
            ρ_anom[i, j, k] = ρij - ρ_mean_k
        end
    end

    # ── Step 1: bp[i,j,k] = ∫ density-anomaly weight (Pa, zero at bottom) ─
    # bp(z) = -∫_bot^z ρ'(z') g dz'  ; integrated with trapezoidal rule on
    # cell-centre-to-cell-centre spacings.
    bp = zeros(Float64, Nx, Ny, Nz)
    @inbounds for j in 1:Ny, i in 1:Nx
        kb = k_bot[i, j]
        kb > Nz && continue
        for k in (kb+1):Nz
            inactive_node(i, j, k, grid, Center(), Center(), Center()) && break
            dza = 0.5 * (Δzᵃᵃᶜ(i, j, k-1, grid) + Δzᵃᵃᶜ(i, j, k, grid))
            bp[i, j, k] = bp[i, j, k-1] -
                          G * 0.5 * (ρ_anom[i, j, k] + ρ_anom[i, j, k-1]) * dza
        end
    end

    # ── Step 2: face-bottom-summed E and 1/H on u- and v-faces ──────────
    # u-face (i+½, j) is between T-cells (i,j) and (ip1,j).
    # v-face (i, j+½) is between T-cells (i,j) and (i,jp1).
    # For each face, k_face = max(k_bot_left, k_bot_right); if that is
    # > Nz the face is entirely below land for at least one side, hence
    # closed (E and invH stay 0).
    E_u    = zeros(Float64, Nx, Ny)        # at (i+½, j)
    E_v    = zeros(Float64, Nx, Ny)        # at (i, j+½)
    invH_u = zeros(Float64, Nx, Ny)
    invH_v = zeros(Float64, Nx, Ny)

    @inbounds for j in 1:Ny, i in 1:Nx
        ip1 = mod1(i + 1, Nx)

        # u-face (i+½, j)
        kf_u = max(k_bot[i, j], k_bot[ip1, j])
        if kf_u <= Nz
            sum_dz = 0.0; sum_bp = 0.0
            for k in kf_u:Nz
                dzk = Δzᵃᵃᶜ(i, j, k, grid)        # Δzᵃᵃᶜ depends on k only, not i,j
                sum_dz += dzk
                sum_bp += 0.5 * (bp[i, j, k] + bp[ip1, j, k]) * dzk
            end
            E_u[i, j]    = sum_bp
            invH_u[i, j] = sum_dz > 0 ? 1.0 / sum_dz : 0.0
        end

        # v-face (i, j+½) — only meaningful if j < Ny (no v-face at Ny+½ inside domain)
        if j < Ny
            kf_v = max(k_bot[i, j], k_bot[i, j+1])
            if kf_v <= Nz
                sum_dz = 0.0; sum_bp = 0.0
                for k in kf_v:Nz
                    dzk = Δzᵃᵃᶜ(i, j, k, grid)
                    sum_dz += dzk
                    sum_bp += 0.5 * (bp[i, j, k] + bp[i, j+1, k]) * dzk
                end
                E_v[i, j]    = sum_bp
                invH_v[i, j] = sum_dz > 0 ? 1.0 / sum_dz : 0.0
            end
        end
    end

    # ── Step 3: J(1/H, E) on T-cells from face-staggered values ─────────
    # ∂(1/H)/∂x and ∂E/∂x at T-cell C use u-faces E and W.
    # ∂(1/H)/∂y and ∂E/∂y at T-cell C use v-faces N and S.
    # Index mapping: u-face (i+½, j) lives at array index (i, j); the W
    # face of T-cell (i, j) is u-face (i-½, j) = array index (im1, j).
    # Similarly v-face (i, j+½) at (i, j); S face of T-cell (i, j) at (i, j-1).
    fill!(ubar_jbar, 0.0)
    @inbounds for j in 2:(Ny-1), i in 1:Nx
        # Skip fully-land cells.
        k_bot[i, j] > Nz && continue

        ip1 = mod1(i + 1, Nx); im1 = mod1(i - 1, Nx)
        dxj = Δxᶜᶜᶜ(i, j, 1, grid)
        dy  = Δyᶜᶜᶜ(i, j, 1, grid)

        d_invH_dx = (invH_u[i,   j  ] - invH_u[im1, j  ]) / dxj
        d_invH_dy = (invH_v[i,   j  ] - invH_v[i,   j-1]) / dy
        dE_dx     = (E_u[i,   j  ]    - E_u[im1, j  ])    / dxj
        dE_dy     = (E_v[i,   j  ]    - E_v[i,   j-1])    / dy

        ubar_jbar[i, j] = d_invH_dx * dE_dy - d_invH_dy * dE_dx
    end
    return ubar_jbar
end

"""
    build_barotropic_solver(grid::AbstractGrid) -> BarotropicSolver

Convenience constructor: extracts ocean mask, bathymetric depth, Coriolis,
and grid spacings from the Oceananigans ImmersedBoundaryGrid, then calls
the array-based assembler.

A column (i,j) is ocean when `bottom_height[i,j] < 0`; H = −bottom_height.
"""
function build_barotropic_solver(grid::AbstractGrid)::BarotropicSolver
    Nx = grid.Nx
    Ny = grid.Ny

    ocean_mask = Matrix{Bool}(undef, Nx, Ny)
    H          = zeros(Float64, Nx, Ny)
    dx         = zeros(Float64, Ny)
    dy_val     = Δyᶜᶜᶜ(1, 1, 1, grid)   # uniform meridional spacing, m

    bh = grid.immersed_boundary.bottom_height   # (Nx, Ny), negative = ocean

    @inbounds for j in 1:Ny
        dx[j] = Δxᶜᶜᶜ(1, j, 1, grid)   # zonal spacing at T-cell centre, m
        for i in 1:Nx
            is_ocean = bh[i, j] < 0.0
            ocean_mask[i, j] = is_ocean
            H[i, j]          = is_ocean ? -bh[i, j] : 0.0
        end
    end

    # Coriolis (with equatorial floor) at T-cell latitudes
    f = [begin
             φ = grid.φᵃᶜᵃ[j]
             fj = 2.0 * (2π / 86164.0) * sind(φ)
             s  = fj >= 0.0 ? 1.0 : -1.0
             s * max(abs(fj), F_MIN)
         end for j in 1:Ny]

    r_bt = compute_rbt(ocean_mask, H, dx, dy_val)
    return build_barotropic_solver(ocean_mask, f, H, dx, dy_val, r_bt)
end

end # module Barotropic
