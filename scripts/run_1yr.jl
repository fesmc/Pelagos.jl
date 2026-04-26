# scripts/run_1yr.jl
#
# Drive the coupled Pelagos OceanModel for 1 model year starting from the
# CLIMBER-X PI steady-state restart, forced by the monthly wind-stress
# climatology in `fake_clim_const_WFDEI-ERA_preind.nc`.  Monthly snapshots of
# T, S, u, v, w, ψ_bt, ρ, and p are written to a single NetCDF file with a
# time axis.
#
# Usage:
#   julia --project=. scripts/run_1yr.jl
#
# Optional environment variables:
#   PELAGOS_RESTART  — path to the CLIMBER-X restart directory
#                       (default: climber-x/restart/pi)
#   PELAGOS_FORCING  — path to the wind-stress NetCDF file
#                       (default: climber-x/input/fake_clim_const_WFDEI-ERA_preind.nc)
#   PELAGOS_OUTPUT   — path to the output NetCDF file
#                       (default: scripts/output/pelagos_1yr_monthly.nc)
#   PELAGOS_DT_DAYS  — timestep in days (default: 1.0)
#   PELAGOS_NO_WIND  — if set to "1", zero out wind stress (diffusion-only baseline)
#
# Stability note: with the CLIMBER-X PI wind-stress climatology, the present
# barotropic ψ solve produces |ψ_bt| ≳ 10⁶ Sv (three orders of magnitude too
# large) and the resulting velocities drive T, S to NaN within the first
# month.  The diffusion-only baseline (PELAGOS_NO_WIND=1) runs cleanly for the
# full year — see project memory.  Treat this script as the harness; the
# physics fix belongs in the barotropic solver, not here.
#
# The "year" used here is the GOLDSTEIN/CLIMBER-X 360-day calendar:
# 12 months × 30 days each.  A snapshot is written at the end of every month.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Pelagos
using Oceananigans
using Oceananigans.Fields: interior
using Oceananigans.Grids: znodes, λnodes, φnodes
using Oceananigans.Operators: Δxᶜᶜᶜ, Δyᶜᶜᶜ
using NCDatasets
using Printf

# ── Configuration ─────────────────────────────────────────────────────────────

const PROJECT_DIR  = abspath(joinpath(@__DIR__, ".."))
const RESTART_DIR  = get(ENV, "PELAGOS_RESTART", joinpath(PROJECT_DIR, "climber-x", "restart", "pi"))
const FORCING_FILE = get(ENV, "PELAGOS_FORCING", joinpath(PROJECT_DIR, "climber-x", "input",
                                                           "fake_clim_const_WFDEI-ERA_preind.nc"))
const OUTPUT_FILE  = get(ENV, "PELAGOS_OUTPUT",  joinpath(@__DIR__, "output", "pelagos_1yr_monthly.nc"))
const DT_DAYS      = parse(Float64, get(ENV, "PELAGOS_DT_DAYS", "1.0"))
const NO_WIND      = get(ENV, "PELAGOS_NO_WIND", "0") == "1"

const DAYS_PER_MONTH = 30
const MONTHS         = 12
const TOTAL_DAYS     = DAYS_PER_MONTH * MONTHS         # 360
const DT_SECONDS     = DT_DAYS * 86400.0
const STEPS_PER_MONTH = round(Int, DAYS_PER_MONTH / DT_DAYS)

isfile(joinpath(RESTART_DIR, "ocn_restart.nc")) ||
    error("No ocn_restart.nc in $RESTART_DIR  (set PELAGOS_RESTART)")
isfile(FORCING_FILE) ||
    error("Forcing file not found: $FORCING_FILE  (set PELAGOS_FORCING)")
mkpath(dirname(OUTPUT_FILE))

@info "Pelagos 1-year run" RESTART_DIR FORCING_FILE OUTPUT_FILE DT_DAYS STEPS_PER_MONTH NO_WIND

# ── Build model and load forcing ──────────────────────────────────────────────

m   = build_ocean_model(RESTART_DIR)
fcg = load_climberx_forcing(FORCING_FILE)   # taux, tauy: (Nx, Ny, 12)

grid = m.grid
Nx, Ny, Nz = grid.Nx, grid.Ny, grid.Nz

size(fcg.taux) == (Nx, Ny, 12) ||
    error("Forcing taux size $(size(fcg.taux)) does not match grid ($Nx, $Ny, 12)")

# ── Coordinate vectors (extracted once) ───────────────────────────────────────

lon       = collect(λnodes(grid, Center(), Center(), Center()))   # °E, T-points
lat       = collect(φnodes(grid, Center(), Center(), Center()))   # °N, T-points
z_centres = collect(znodes(grid, Center(), Center(), Center()))   # m, T-points (negative)
z_faces   = collect(znodes(grid, Center(), Center(), Face()))     # m, w-points (negative)

# ── Diagnostic helpers ────────────────────────────────────────────────────────

# Bathymetric geometry needed for the diagnostic ψ_bt solve.
bh = grid.immersed_boundary.bottom_height        # (Nx, Ny), negative = ocean
ocean_mask_2d = Matrix{Bool}(undef, Nx, Ny)
H_2d          = zeros(Float64, Nx, Ny)
for j in 1:Ny, i in 1:Nx
    is_ocean = bh[i, j] < 0.0
    ocean_mask_2d[i, j] = is_ocean
    H_2d[i, j]          = is_ocean ? -bh[i, j] : 0.0
end

dx_vec = [Δxᶜᶜᶜ(1, j, 1, grid) for j in 1:Ny]
dy_val = Δyᶜᶜᶜ(1, 1, 1, grid)

const F_FLOOR = 5e-6
f_vec = [begin
             φ  = grid.φᵃᶜᵃ[j]
             fj = 2.0 * (2π / 86164.0) * sind(φ)
             s  = fj >= 0.0 ? 1.0 : -1.0
             s * max(abs(fj), F_FLOOR)
         end for j in 1:Ny]

snapshot_psi_bt(model, τx, τy) =
    solve_barotropic!(model.bar_solver, τx, τy, H_2d, f_vec, dx_vec, dy_val, ocean_mask_2d)

function snapshot_density(T, S, p)
    ρ = zeros(Float64, size(T))
    @inbounds for k in axes(T,3), j in axes(T,2), i in axes(T,1)
        Sij = max(S[i,j,k], 0.0)
        ρ[i,j,k] = seawater_density(T[i,j,k], Sij, max(p[i,j,k], 0.0))
    end
    return ρ
end

# Velocities live on staggered C-grid faces.  Interpolate to T-points so output
# has uniform (Nx, Ny, Nz) shape.  Periodic in λ; bounded in φ and z.
function u_to_tpoint(u_face::AbstractArray{<:Real,3})
    Nx, Ny, Nz = size(u_face)
    out = zeros(Float64, Nx, Ny, Nz)
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        ip1       = mod1(i + 1, Nx)
        out[i,j,k] = 0.5 * (u_face[i, j, k] + u_face[ip1, j, k])
    end
    return out
end

function v_to_tpoint(v_face::AbstractArray{<:Real,3})
    _, Ny1, Nz = size(v_face)
    Ny  = Ny1 - 1
    Nx  = size(v_face, 1)
    out = zeros(Float64, Nx, Ny, Nz)
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        out[i,j,k] = 0.5 * (v_face[i, j, k] + v_face[i, j+1, k])
    end
    return out
end

function w_to_tpoint(w_face::AbstractArray{<:Real,3})
    Nx, Ny, Nz1 = size(w_face)
    Nz  = Nz1 - 1
    out = zeros(Float64, Nx, Ny, Nz)
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        out[i,j,k] = 0.5 * (w_face[i, j, k] + w_face[i, j, k+1])
    end
    return out
end

# ── NetCDF output setup ───────────────────────────────────────────────────────

isfile(OUTPUT_FILE) && rm(OUTPUT_FILE)

ds = NCDatasets.Dataset(OUTPUT_FILE, "c")
NCDatasets.defDim(ds, "lon",  Nx)
NCDatasets.defDim(ds, "lat",  Ny)
NCDatasets.defDim(ds, "z",    Nz)
NCDatasets.defDim(ds, "zf",   Nz + 1)
NCDatasets.defDim(ds, "time", MONTHS)

ds.attrib["title"]       = "Pelagos.jl 1-year monthly diagnostic run"
ds.attrib["calendar"]    = "360_day (12 × 30-day months)"
ds.attrib["dt_days"]     = DT_DAYS
ds.attrib["restart_dir"] = RESTART_DIR
ds.attrib["forcing_file"] = FORCING_FILE

function defvar!(ds, name, T, dims, units, long_name)
    v = NCDatasets.defVar(ds, name, T, dims)
    v.attrib["units"]     = units
    v.attrib["long_name"] = long_name
    return v
end

v_lon  = defvar!(ds, "lon",  Float32, ("lon",),  "degrees_east",  "longitude")
v_lat  = defvar!(ds, "lat",  Float32, ("lat",),  "degrees_north", "latitude")
v_z    = defvar!(ds, "z",    Float32, ("z",),    "m",             "T-point depth (negative downward)")
v_zf   = defvar!(ds, "zf",   Float32, ("zf",),   "m",             "w-point interface depth (negative downward)")
v_time = defvar!(ds, "time", Float32, ("time",), "days since simulation start", "month-end time")

v_T   = defvar!(ds, "T",      Float32, ("lon","lat","z","time"),  "degC",  "potential temperature")
v_S   = defvar!(ds, "S",      Float32, ("lon","lat","z","time"),  "psu",   "salinity")
v_u   = defvar!(ds, "u",      Float32, ("lon","lat","z","time"),  "m s-1", "zonal velocity (C-grid east face)")
v_v   = defvar!(ds, "v",      Float32, ("lon","lat","z","time"),  "m s-1", "meridional velocity (C-grid north face)")
v_w   = defvar!(ds, "w",      Float32, ("lon","lat","z","time"),  "m s-1", "vertical velocity (T-point, averaged from z-faces)")
v_w_f = defvar!(ds, "w_face", Float32, ("lon","lat","zf","time"), "m s-1", "vertical velocity at z-faces")
v_ψ   = defvar!(ds, "psi_bt", Float32, ("lon","lat","time"),      "m3 s-1","barotropic streamfunction (volume transport)")
v_ρ   = defvar!(ds, "rho",    Float32, ("lon","lat","z","time"),  "kg m-3","UNESCO in-situ density")
v_p   = defvar!(ds, "p",      Float32, ("lon","lat","z","time"),  "bar",   "hydrostatic pressure")

v_mask = defvar!(ds, "ocean_mask",    Int8,   ("lon","lat"), "1", "1 = ocean column, 0 = land")
v_H    = defvar!(ds, "bathy_depth",   Float32,("lon","lat"), "m", "ocean column depth (positive)")
v_τx   = defvar!(ds, "tau_x_monthly", Float32,("lon","lat","time"), "N m-2",
                                              "zonal wind stress applied during month")
v_τy   = defvar!(ds, "tau_y_monthly", Float32,("lon","lat","time"), "N m-2",
                                              "meridional wind stress applied during month")

v_lon[:]  = Float32.(lon)
v_lat[:]  = Float32.(lat)
v_z[:]    = Float32.(z_centres)
v_zf[:]   = Float32.(z_faces)
v_mask[:] = Int8.(ocean_mask_2d)
v_H[:]    = Float32.(H_2d)

# ── Time loop ─────────────────────────────────────────────────────────────────

T_field = m.tracer_model.tracers.T
S_field = m.tracer_model.tracers.S
u_field = m.tracer_model.velocities.u
v_field = m.tracer_model.velocities.v
w_field = m.tracer_model.velocities.w

function run_year!(m, fcg, ocean_mask_2d, Nz,
                   T_field, S_field, u_field, v_field, w_field,
                   v_T, v_S, v_u, v_v, v_w, v_w_f, v_p, v_ρ, v_ψ, v_τx, v_τy, v_time)
    day = 0.0
    t0  = time()
    @info "Starting integration" total_days=TOTAL_DAYS dt_days=DT_DAYS

    for month in 1:MONTHS
        τx = NO_WIND ? zeros(size(fcg.taux, 1), size(fcg.taux, 2)) : fcg.taux[:, :, month]
        τy = NO_WIND ? zeros(size(fcg.tauy, 1), size(fcg.tauy, 2)) : fcg.tauy[:, :, month]

        for _ in 1:STEPS_PER_MONTH
            step_ocean!(m, τx, τy, DT_SECONDS)
            day += DT_DAYS
        end

        Tarr   = Array(interior(T_field))
        Sarr   = Array(interior(S_field))
        u_face = Array(interior(u_field))
        v_face = Array(interior(v_field))
        w_face = Array(interior(w_field))
        parr   = Array(interior(m.p))
        u_T    = u_to_tpoint(u_face)
        v_T_   = v_to_tpoint(v_face)
        w_T    = w_to_tpoint(w_face)
        ψ      = snapshot_psi_bt(m, τx, τy)
        ρ      = snapshot_density(Tarr, Sarr, parr)

        v_T[:,:,:,month]   = Float32.(Tarr)
        v_S[:,:,:,month]   = Float32.(Sarr)
        v_u[:,:,:,month]   = Float32.(u_T)
        v_v[:,:,:,month]   = Float32.(v_T_)
        v_w[:,:,:,month]   = Float32.(w_T)
        v_w_f[:,:,:,month] = Float32.(w_face)
        v_p[:,:,:,month]   = Float32.(parr)
        v_ρ[:,:,:,month]   = Float32.(ρ)
        v_ψ[:,:,month]     = Float32.(ψ)
        v_τx[:,:,month]    = Float32.(τx)
        v_τy[:,:,month]    = Float32.(τy)
        v_time[month]      = Float32(day)

        om3     = repeat(ocean_mask_2d, 1, 1, Nz)
        T_ocean = Tarr[om3]
        S_ocean = Sarr[om3]
        @info @sprintf("month %2d  day %3d  T̄=%.3f  T-range=[%.2f,%.2f]  S̄=%.3f  S-range=[%.2f,%.2f]  ψ_bt±=[%.1f,%.1f] Sv  elapsed=%.1fs",
            month, Int(day),
            sum(T_ocean)/length(T_ocean), minimum(T_ocean), maximum(T_ocean),
            sum(S_ocean)/length(S_ocean), minimum(S_ocean), maximum(S_ocean),
            minimum(ψ)/1e6, maximum(ψ)/1e6,
            time() - t0)
    end
    return time() - t0
end

wallclock = run_year!(m, fcg, ocean_mask_2d, Nz,
                      T_field, S_field, u_field, v_field, w_field,
                      v_T, v_S, v_u, v_v, v_w, v_w_f, v_p, v_ρ, v_ψ, v_τx, v_τy, v_time)

close(ds)
@info "Done" output=OUTPUT_FILE wallclock=wallclock
