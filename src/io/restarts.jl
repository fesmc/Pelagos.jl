# Restart I/O for Pelagos.jl.
#
# Two distinct restart formats are handled:
#
#   Pelagos native  — write_restart / read_restart
#     Minimal NetCDF checkpoints written by Pelagos itself (T, S, time).
#
#   CLIMBER-X       — load_climberx_restart / load_climberx_forcing
#     Read the CLIMBER-X ocn_restart.nc to initialise Pelagos from a
#     CLIMBER-X steady-state.  Also reads monthly wind stress climatology
#     from the atmospheric forcing file used by the default PI control run.
#
# CLIMBER-X grid convention
# ─────────────────────────
# Tracers (T, S) are on T-points: (nlon=72, nlat=36, nz=23).
# Velocities are on B-grid corners: u(3, nlon1=73, nlat1=37, nz=23) where
# dim-1 carries (u, v, w).  Converting to the Oceananigans C-grid requires
# averaging adjacent corner values:
#   u_C[i, j, k] = (u_B[i, j, k] + u_B[i, j+1, k]) / 2   (average in φ)
#   v_C[i, j, k] = (v_B[i, j, k] + v_B[i+1, j, k]) / 2   (average in λ, periodic)
#
# Layer ordering: index 1 = deepest layer (~5000 m), index 23 = shallowest (~5 m).
# This matches Oceananigans' bottom-to-top convention for k.

module Restarts

using NCDatasets

export write_restart, read_restart
export load_climberx_restart, load_climberx_forcing
export bgrid_u_to_cgrid, bgrid_v_to_cgrid

# ── Pelagos native restart ─────────────────────────────────────────────────────

"""
    write_restart(filename, T, S, time_yr)

Write a Pelagos restart file with tracer fields and current model time.
"""
function write_restart(filename::String,
                       T       ::AbstractArray{Float64,3},
                       S       ::AbstractArray{Float64,3},
                       time_yr ::Float64)
    nlon, nlat, nz = size(T)
    NCDatasets.Dataset(filename, "c") do ds
        NCDatasets.defDim(ds, "nlon", nlon)
        NCDatasets.defDim(ds, "nlat", nlat)
        NCDatasets.defDim(ds, "nz",   nz)

        v_time = NCDatasets.defVar(ds, "time", Float64, ())
        v_time.attrib["units"] = "years"
        v_time[] = time_yr

        v_T = NCDatasets.defVar(ds, "T", Float64, ("nlon", "nlat", "nz"))
        v_T.attrib["units"] = "degC"
        v_T[:] = T

        v_S = NCDatasets.defVar(ds, "S", Float64, ("nlon", "nlat", "nz"))
        v_S.attrib["units"] = "psu"
        v_S[:] = S
    end
    return filename
end

"""
    read_restart(filename) -> (T, S, time_yr)

Read a Pelagos native restart file.
"""
function read_restart(filename::String)
    NCDatasets.Dataset(filename, "r") do ds
        T       = Float64.(ds["T"][:])
        S       = Float64.(ds["S"][:])
        time_yr = Float64(ds["time"][])
        return (T, S, time_yr)
    end
end

# ── CLIMBER-X restart reader ───────────────────────────────────────────────────

"""
    load_climberx_restart(restart_dir) -> NamedTuple

Load the CLIMBER-X ocean steady-state from `<restart_dir>/ocn_restart.nc`.

Returned fields (all Float64, bottom-to-top layer ordering):
- `T`     : temperature (72, 36, 23), °C
- `S`     : salinity    (72, 36, 23), psu
- `u_B`   : zonal velocity on B-grid corners     (73, 37, 23), m s⁻¹
- `v_B`   : meridional velocity on B-grid corners (73, 37, 23), m s⁻¹
- `w_B`   : vertical velocity on B-grid corners   (73, 37, 23), m s⁻¹
- `f_ocn` : ocean fraction mask (72, 36)
- `k1_pot`: bottom layer index  (72, 36)
- `zw`    : layer interfaces (24,), m, negative, bottom→top
"""
function load_climberx_restart(restart_dir::String)
    file = joinpath(restart_dir, "ocn_restart.nc")
    isfile(file) || error("ocn_restart.nc not found in \"$restart_dir\"")

    NCDatasets.Dataset(file, "r") do ds
        # ts is stored flattened (lon*lat*nz*ntrc), reshape and extract T/S
        ts = reshape(Float64.(ds["ts"][:]), 72, 36, 23, 4)
        T  = ts[:, :, :, 1]
        S  = ts[:, :, :, 2]

        # u is stored flattened (3*lon1*lat1*nz), dim 1 = (u, v, w)
        u_all = reshape(Float64.(ds["u"][:]), 3, 73, 37, 23)
        u_B   = u_all[1, :, :, :]    # (73, 37, 23)
        v_B   = u_all[2, :, :, :]
        w_B   = u_all[3, :, :, :]

        f_ocn = reshape(Float64.(ds["f_ocn"][:]), 72, 36)
        k1    = reshape(Int.(ds["k1_pot"][:]),    72, 36)
        zw    = Float64.(ds["zw"][:])             # (24,)

        (T=T, S=S, u_B=u_B, v_B=v_B, w_B=w_B, f_ocn=f_ocn, k1_pot=k1, zw=zw)
    end
end

# ── CLIMBER-X forcing reader ───────────────────────────────────────────────────

"""
    load_climberx_forcing(forcing_file) -> NamedTuple

Load monthly wind stress climatology from a CLIMBER-X atmospheric forcing file
(e.g., `climber-x/input/fake_clim_const_WFDEI-ERA_preind.nc`).

Returned fields (Float64):
- `taux` : (72, 36, 12) zonal wind stress,      N m⁻²
- `tauy` : (72, 36, 12) meridional wind stress,  N m⁻²
"""
function load_climberx_forcing(forcing_file::String)
    isfile(forcing_file) || error("forcing file not found: \"$forcing_file\"")
    NCDatasets.Dataset(forcing_file, "r") do ds
        taux = Float64.(ds["taux"][:, :, :])   # (72, 36, 12)
        tauy = Float64.(ds["tauy"][:, :, :])
        (taux=taux, tauy=tauy)
    end
end

# ── B-grid → C-grid velocity interpolation ────────────────────────────────────

"""
    bgrid_u_to_cgrid(u_B) -> u_C

Convert zonal velocity from CLIMBER-X B-grid corners (73, 37, nz) to
Oceananigans C-grid east faces (73, 36, nz) by averaging adjacent φ values.
"""
function bgrid_u_to_cgrid(u_B::AbstractArray{Float64,3})
    _, nlat1, nz = size(u_B)   # (73, 37, nz)
    nlat = nlat1 - 1           # 36
    u_C  = similar(u_B, 73, nlat, nz)
    @inbounds for k in 1:nz, j in 1:nlat
        u_C[:, j, k] = (u_B[:, j, k] .+ u_B[:, j+1, k]) ./ 2
    end
    return u_C
end

"""
    bgrid_v_to_cgrid(v_B) -> v_C

Convert meridional velocity from CLIMBER-X B-grid corners (73, 37, nz) to
Oceananigans C-grid north faces (72, 37, nz) by averaging adjacent λ values
(periodic in longitude).
"""
function bgrid_v_to_cgrid(v_B::AbstractArray{Float64,3})
    nlon1, _, nz = size(v_B)   # (73, 37, nz)
    nlon  = nlon1 - 1          # 72
    v_C   = similar(v_B, nlon, 37, nz)
    @inbounds for k in 1:nz, i in 1:nlon
        i_next = mod1(i + 1, nlon1)
        v_C[i, :, k] = (v_B[i, :, k] .+ v_B[i_next, :, k]) ./ 2
    end
    return v_C
end

end # module Restarts
