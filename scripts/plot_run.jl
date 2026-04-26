# scripts/plot_run.jl
#
# Visualize a Pelagos.jl 1-year monthly NetCDF output (as produced by
# scripts/run_1yr.jl) using CairoMakie.
#
# Produces a set of PNG figures in <input_dir>/figures/:
#   fig0_summary.png         — one-page overview (SST, ψ_bt, T section, time series)
#   fig1_surface_maps.png    — annual-mean SST, SSS, ψ_bt, surface speed + arrows
#   fig2_zonal_sections.png  — zonal-mean T, S, ρ versus depth × latitude
#   fig3_zonal_profiles.png  — zonal-mean SST, SSS, ψ_bt versus latitude
#   fig4_moc.png             — meridional overturning streamfunction
#   fig5_timeseries.png      — global-mean T, S, peak overturning at 26°N
#
# Usage:
#   julia --project=scripts scripts/plot_run.jl [INPUT.nc] [OUTDIR]
#
# Defaults:
#   INPUT  = scripts/output/pelagos_1yr_monthly.nc
#   OUTDIR = <dirname of INPUT>/figures
#
# The first run will instantiate the script-local Project.toml (CairoMakie +
# NCDatasets + ColorSchemes); subsequent runs reuse the precompiled cache.

using Pkg
Pkg.activate(@__DIR__)
try
    using CairoMakie, NCDatasets, ColorSchemes, Statistics
catch
    Pkg.instantiate()
    using CairoMakie, NCDatasets, ColorSchemes, Statistics
end

const Sv = 1e6        # 1 Sverdrup = 10⁶ m³ s⁻¹
const R_EARTH = 6.371e6

# ── Argument parsing ─────────────────────────────────────────────────────────

const DEFAULT_INPUT = joinpath(@__DIR__, "output", "pelagos_1yr_monthly.nc")
input_path  = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_INPUT
isfile(input_path) || error("input file not found: $input_path")
output_dir  = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(input_path), "figures")
mkpath(output_dir)
@info "Pelagos plot_run" input_path output_dir

# ── Load data ────────────────────────────────────────────────────────────────

ds = NCDatasets.Dataset(input_path, "r")

lon       = Float64.(ds["lon"][:])
lat       = Float64.(ds["lat"][:])
z_centres = Float64.(ds["z"][:])
z_faces   = Float64.(ds["zf"][:])
times     = Float64.(ds["time"][:])

T   = Float64.(ds["T"][:, :, :, :])        # (lon, lat, z, time)
S   = Float64.(ds["S"][:, :, :, :])
u   = Float64.(ds["u"][:, :, :, :])
v   = Float64.(ds["v"][:, :, :, :])
w   = Float64.(ds["w"][:, :, :, :])
ρ   = Float64.(ds["rho"][:, :, :, :])
ψbt = Float64.(ds["psi_bt"][:, :, :])      # (lon, lat, time)
ocean_mask = ds["ocean_mask"][:, :] .!= 0
H_depth    = Float64.(ds["bathy_depth"][:, :])

close(ds)

Nx, Ny, Nz = length(lon), length(lat), length(z_centres)
Ntime      = length(times)
@info "Grid" Nx Ny Nz Ntime

# Apply land mask everywhere
mask3 = repeat(ocean_mask, 1, 1, Nz)
mask4 = repeat(mask3, 1, 1, 1, Ntime)
mask_psi4 = repeat(ocean_mask, 1, 1, Ntime)

T[.!mask4]   .= NaN
S[.!mask4]   .= NaN
u[.!mask4]   .= NaN
v[.!mask4]   .= NaN
w[.!mask4]   .= NaN
ρ[.!mask4]   .= NaN
ψbt[.!mask_psi4] .= NaN

# ── Helpers ──────────────────────────────────────────────────────────────────

annual_mean(A; dims) = dropdims(mean(skipmissing_array(A); dims=dims); dims=dims)

# Mean ignoring NaNs along given dim, returns array reduced over `dims`.
function nanmean(A; dims)
    s = sum(x -> isnan(x) ? 0.0 : x, A; dims=dims)
    n = sum(x -> isnan(x) ? 0   : 1, A; dims=dims)
    out = s ./ n
    out[n .== 0] .= NaN
    return out
end
nanmean1(A; dims) = dropdims(nanmean(A; dims=dims); dims=dims)

skipmissing_array(A) = A   # placeholder; we do explicit nanmean above

# Cell area (m²) on the lat-lon grid
function cell_area(lon, lat)
    dlon = deg2rad(abs(lon[2] - lon[1]))
    dlat = deg2rad(abs(lat[2] - lat[1]))
    cosφ = cos.(deg2rad.(lat))
    return [R_EARTH^2 * dlon * dlat * cosφ[j] for i in eachindex(lon), j in eachindex(lat)]
end

const AREA   = cell_area(lon, lat)
const DZ     = diff(z_faces)              # (Nz,)
const DX_LAT = R_EARTH * deg2rad(abs(lon[2]-lon[1])) .* cos.(deg2rad.(lat))   # (Ny,)

# Build a divergent range centred at zero — never returns a degenerate (a, a)
function div_norm(field)
    finite = filter(isfinite, field)
    isempty(finite) && return (-1.0, 1.0)
    a = maximum(abs, finite)
    a == 0 && return (-1.0, 1.0)
    return (-a, a)
end

# Robust extrema — guards against all-NaN or constant fields (would otherwise
# crash CairoMakie's continuous color sampling).
function safe_range(field)
    finite = filter(isfinite, field)
    isempty(finite) && return (0.0, 1.0)
    lo, hi = extrema(finite)
    if lo == hi
        eps = abs(lo) > 0 ? 0.01 * abs(lo) : 1.0
        return (lo - eps, hi + eps)
    end
    return (lo, hi)
end

# ── Figure 1: surface maps (annual mean) ─────────────────────────────────────

function fig_surface_maps()
    sst    = nanmean1(T[:, :, end:end, :]; dims=4)[:, :, 1]   # (Nx, Ny)
    sss    = nanmean1(S[:, :, end:end, :]; dims=4)[:, :, 1]
    ψam    = nanmean1(ψbt; dims=3) ./ Sv
    usurf  = nanmean1(u[:, :, end:end, :]; dims=4)[:, :, 1]
    vsurf  = nanmean1(v[:, :, end:end, :]; dims=4)[:, :, 1]
    speed  = sqrt.(usurf.^2 .+ vsurf.^2)

    fig = Figure(size = (1200, 950))

    function add_map!(ax, field; title, cmap, divergent=false, label="")
        ax.title = title
        ax.xlabel = "Longitude (°E)"
        ax.ylabel = "Latitude (°N)"
        if divergent
            lo, hi = div_norm(field)
            hm = heatmap!(ax, lon, lat, field; colormap=cmap, colorrange=(lo, hi))
        else
            hm = heatmap!(ax, lon, lat, field; colormap=cmap)
        end
        Colorbar(fig[:, 0], hm; label=label, vertical=false)   # placeholder
        return hm
    end

    ax1 = Axis(fig[1, 1]; title="Annual mean SST",
               xlabel="°E", ylabel="°N")
    hm1 = heatmap!(ax1, lon, lat, sst; colormap=:RdYlBu_4, colorrange=safe_range(sst))
    Colorbar(fig[1, 2], hm1, label="°C")

    ax2 = Axis(fig[1, 3]; title="Annual mean SSS", xlabel="°E", ylabel="°N")
    hm2 = heatmap!(ax2, lon, lat, sss; colormap=:viridis, colorrange=safe_range(sss))
    Colorbar(fig[1, 4], hm2, label="psu")

    ax3 = Axis(fig[2, 1]; title="Barotropic ψ (annual mean)", xlabel="°E", ylabel="°N")
    lo, hi = div_norm(ψam)
    hm3 = heatmap!(ax3, lon, lat, ψam; colormap=:RdBu, colorrange=(lo, hi))
    Colorbar(fig[2, 2], hm3, label="Sv")

    ax4 = Axis(fig[2, 3]; title="Surface speed + velocity (annual mean)",
               xlabel="°E", ylabel="°N")
    hm4 = heatmap!(ax4, lon, lat, speed; colormap=:magma,
                   colorrange=safe_range(speed))
    # Subsampled velocity arrows (every 2 cells)
    is = 1:2:Nx; js = 1:2:Ny
    pts = [Point2f(lon[i], lat[j]) for i in is, j in js if isfinite(usurf[i,j])]
    dirs= [Vec2f(usurf[i,j], vsurf[i,j]) for i in is, j in js if isfinite(usurf[i,j])]
    if !isempty(pts)
        arrows2d!(ax4, pts, dirs; lengthscale = 4.0, color=:white)
    end
    Colorbar(fig[2, 4], hm4, label="m/s")

    Label(fig[0, 1:4], "Surface diagnostics — annual mean";
          fontsize=18, padding=(0,0,0,10))
    out = joinpath(output_dir, "fig1_surface_maps.png")
    save(out, fig); @info "wrote" out
end

# ── Figure 2: zonal-mean depth-latitude sections ─────────────────────────────

function fig_zonal_sections()
    Tann = nanmean1(T; dims=4)               # (Nx, Ny, Nz)
    Sann = nanmean1(S; dims=4)
    rann = nanmean1(ρ; dims=4)
    Tzm  = nanmean1(Tann; dims=1)            # (Ny, Nz)
    Szm  = nanmean1(Sann; dims=1)
    rzm  = nanmean1(rann; dims=1)

    fig = Figure(size = (1500, 460))

    function panel!(col, M, title, cmap, label)
        ax = Axis(fig[1, col]; title, xlabel="Latitude (°N)",
                  ylabel = col == 1 ? "Depth (m)" : "")
        hm = heatmap!(ax, lat, z_centres, M;
                      colormap=cmap,
                      colorrange = safe_range(M))
        Colorbar(fig[2, col], hm; label=label, vertical=false, flipaxis=false)
    end

    panel!(1, Tzm, "Zonal-mean T",     :RdYlBu_4, "°C")
    panel!(2, Szm, "Zonal-mean S",     :viridis,  "psu")
    panel!(3, rzm, "Zonal-mean ρ",     :cividis,  "kg m⁻³")

    Label(fig[0, 1:3], "Zonal-mean depth–latitude sections (annual mean)";
          fontsize=18, padding=(0,0,0,10))
    out = joinpath(output_dir, "fig2_zonal_sections.png")
    save(out, fig); @info "wrote" out
end

# ── Figure 3: zonal-mean profiles vs latitude ────────────────────────────────

function fig_zonal_profiles()
    sst = nanmean1(T[:, :, end:end, :]; dims=4)[:, :, 1]   # (Nx, Ny)
    sss = nanmean1(S[:, :, end:end, :]; dims=4)[:, :, 1]
    ψam = nanmean1(ψbt; dims=3) ./ Sv

    sst_lat = nanmean1(sst; dims=1)
    sss_lat = nanmean1(sss; dims=1)
    ψ_lat   = nanmean1(ψam; dims=1)

    fig = Figure(size = (1400, 380))
    ax1 = Axis(fig[1,1]; title="Zonal-mean SST", xlabel="Latitude (°N)", ylabel="°C")
    lines!(ax1, lat, sst_lat; color=:firebrick, linewidth=2)
    ax2 = Axis(fig[1,2]; title="Zonal-mean SSS", xlabel="Latitude (°N)", ylabel="psu")
    lines!(ax2, lat, sss_lat; color=:navy, linewidth=2)
    ax3 = Axis(fig[1,3]; title="Zonal-mean ψ_bt", xlabel="Latitude (°N)", ylabel="Sv")
    lines!(ax3, lat, ψ_lat; color=:purple, linewidth=2)
    for ax in (ax1, ax2, ax3); ax.xgridvisible[]=true; ax.ygridvisible[]=true; end
    Label(fig[0, 1:3], "Zonal-mean profiles vs latitude (annual mean)";
          fontsize=18, padding=(0,0,0,10))
    out = joinpath(output_dir, "fig3_zonal_profiles.png")
    save(out, fig); @info "wrote" out
end

# ── Figure 4: meridional overturning streamfunction ─────────────────────────

"""Compute meridional overturning streamfunction ψ_moc(lat, z) in Sv from
zonally integrated v at T-points and a vertical cumulative integral upward."""
function compute_moc(v_array)
    # v_array: (Nx, Ny, Nz) annual mean — NaNs over land
    v_clean = replace(v_array, NaN => 0.0)
    v_zonal = zeros(Float64, Ny, Nz)
    for k in 1:Nz, j in 1:Ny
        # zonal integral with metric DX_LAT[j] (uniform Δλ at this latitude)
        s = 0.0
        @inbounds for i in 1:Nx
            s += v_clean[i, j, k]
        end
        v_zonal[j, k] = s * DX_LAT[j]
    end
    # cumulative integral upward from the bottom (k=1 = deepest in our convention)
    ψ = zeros(Float64, Ny, Nz)
    for j in 1:Ny
        acc = 0.0
        for k in 1:Nz
            acc += v_zonal[j, k] * DZ[k]
            ψ[j, k] = acc
        end
    end
    return ψ ./ Sv     # → Sv
end

function fig_overturning()
    vann = nanmean1(v; dims=4)
    ψmoc = compute_moc(vann)

    fig = Figure(size = (900, 500))
    ax = Axis(fig[1,1]; title="Meridional overturning streamfunction (annual mean)",
              xlabel="Latitude (°N)", ylabel="Depth (m)")
    lo, hi = div_norm(ψmoc)
    hm = heatmap!(ax, lat, z_centres, ψmoc; colormap=:RdBu, colorrange=(lo, hi))
    contour!(ax, lat, z_centres, ψmoc; levels=10, color=(:black, 0.4), linewidth=0.5)
    Colorbar(fig[1, 2], hm; label="Sv")
    out = joinpath(output_dir, "fig4_moc.png")
    save(out, fig); @info "wrote" out
    return ψmoc
end

# ── Figure 5: time series ────────────────────────────────────────────────────

function fig_timeseries()
    # Volume-weighted global means (use only ocean cells, weight by AREA*DZ)
    weights3 = [AREA[i,j] * DZ[k] for i in 1:Nx, j in 1:Ny, k in 1:Nz]
    weights3 = weights3 .* mask3              # zero out land

    Tg = zeros(Ntime)
    Sg = zeros(Ntime)
    for t in 1:Ntime
        Tt = T[:, :, :, t]
        St = S[:, :, :, t]
        Tt[.!mask3] .= 0.0
        St[.!mask3] .= 0.0
        Tg[t] = sum(Tt .* weights3) / sum(weights3)
        Sg[t] = sum(St .* weights3) / sum(weights3)
    end

    # AMOC peak at 26°N — recompute MOC for each month
    j26 = argmin(abs.(lat .- 26.0))
    amoc26 = zeros(Ntime)
    for t in 1:Ntime
        ψmoc = compute_moc(v[:, :, :, t])
        col  = ψmoc[j26, :]
        finite = filter(isfinite, col)
        amoc26[t] = isempty(finite) ? NaN : maximum(finite)
    end

    fig = Figure(size = (1400, 380))
    ax1 = Axis(fig[1,1]; title="Global-mean T",   xlabel="Days from start", ylabel="°C")
    scatterlines!(ax1, times, Tg)
    ax2 = Axis(fig[1,2]; title="Global-mean S",   xlabel="Days from start", ylabel="psu")
    scatterlines!(ax2, times, Sg; color=:navy)
    ax3 = Axis(fig[1,3]; title="Peak overturning at $(round(Int,lat[j26]))°N",
               xlabel="Days from start", ylabel="Sv")
    scatterlines!(ax3, times, amoc26; color=:purple)
    Label(fig[0, 1:3], "Monthly evolution"; fontsize=18, padding=(0,0,0,10))
    out = joinpath(output_dir, "fig5_timeseries.png")
    save(out, fig); @info "wrote" out
end

# ── Figure 0: one-page summary ───────────────────────────────────────────────

function fig_summary()
    sst = nanmean1(T[:, :, end:end, :]; dims=4)[:, :, 1]
    ψam = nanmean1(ψbt; dims=3) ./ Sv
    Tann = nanmean1(T; dims=4)
    Tzm  = nanmean1(Tann; dims=1)

    fig = Figure(size = (1300, 900))

    ax1 = Axis(fig[1,1]; title="Annual-mean SST",   xlabel="°E", ylabel="°N")
    hm1 = heatmap!(ax1, lon, lat, sst; colormap=:RdYlBu_4,
                   colorrange = safe_range(sst))
    Colorbar(fig[1,2], hm1, label="°C")

    ax2 = Axis(fig[1,3]; title="Barotropic ψ (annual mean)", xlabel="°E", ylabel="°N")
    lo, hi = div_norm(ψam)
    hm2 = heatmap!(ax2, lon, lat, ψam; colormap=:RdBu, colorrange=(lo, hi))
    Colorbar(fig[1,4], hm2, label="Sv")

    ax3 = Axis(fig[2,1:2]; title="Zonal-mean T (annual)",
               xlabel="Latitude (°N)", ylabel="Depth (m)")
    hm3 = heatmap!(ax3, lat, z_centres, Tzm; colormap=:RdYlBu_4,
                   colorrange=safe_range(Tzm))
    Colorbar(fig[2,3], hm3, label="°C")

    # surface time series
    weights2 = AREA .* ocean_mask
    sst_ts = [sum(replace(T[:, :, end, t] .* weights2, NaN=>0.0)) / sum(weights2) for t in 1:Ntime]
    sss_ts = [sum(replace(S[:, :, end, t] .* weights2, NaN=>0.0)) / sum(weights2) for t in 1:Ntime]
    ax4 = Axis(fig[2,4]; title="Surface time series",
               xlabel="Days from start", ylabel="SST (°C)")
    scatterlines!(ax4, times, sst_ts; color=:firebrick, label="SST")
    ax4r = Axis(fig[2,4]; ylabel="SSS (psu)", yaxisposition=:right,
                xticklabelsvisible=false, xticksvisible=false,
                xgridvisible=false, ygridvisible=false)
    scatterlines!(ax4r, times, sss_ts; color=:navy, marker=:rect, label="SSS")
    hidexdecorations!(ax4r); ax4r.xlabel=""

    Label(fig[0, 1:4], "Pelagos.jl 1-year run summary"; fontsize=20,
          padding=(0,0,0,10))
    out = joinpath(output_dir, "fig0_summary.png")
    save(out, fig); @info "wrote" out
end

# ── Run all ──────────────────────────────────────────────────────────────────

fig_summary()
fig_surface_maps()
fig_zonal_sections()
fig_zonal_profiles()
fig_overturning()
fig_timeseries()

@info "Done"
