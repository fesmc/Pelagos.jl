# NetCDF output for Pelagos.jl.
#
# Writes ocean state (T, S, u, v, w, psi_bt) to NetCDF files compatible with
# the CLIMBER-X output conventions so that validation scripts can diff directly
# against reference fixtures.

module Output

export write_snapshot

"""
    write_snapshot(filename, state, grid_info, time_yr)

Write an ocean state snapshot to a NetCDF file.

# Arguments
- `filename`  : output file path (.nc)
- `state`     : NamedTuple with fields: T, S, u, v, w, psi_bt, rho, pressure
- `grid_info` : NamedTuple with: lon, lat, z_centres, z_faces
- `time_yr`   : model time in years (written as a coordinate)
"""
function write_snapshot(filename ::String,
                        state    ,
                        grid_info,
                        time_yr  ::Float64)
    NCDatasets = _load_ncdatasets()

    NCDatasets.Dataset(filename, "c") do ds
        nlon = length(grid_info.lon)
        nlat = length(grid_info.lat)
        nz   = length(grid_info.z_centres)

        # Dimensions
        NCDatasets.defDim(ds, "lon",  nlon)
        NCDatasets.defDim(ds, "lat",  nlat)
        NCDatasets.defDim(ds, "z",    nz)
        NCDatasets.defDim(ds, "zf",   nz+1)
        NCDatasets.defDim(ds, "time", 1)

        # Coordinate variables
        _defvar(ds, "lon",  Float32, ("lon",),  grid_info.lon;  units="degrees_east")
        _defvar(ds, "lat",  Float32, ("lat",),  grid_info.lat;  units="degrees_north")
        _defvar(ds, "z",    Float32, ("z",),    grid_info.z_centres; units="m")
        _defvar(ds, "time", Float32, ("time",), [Float32(time_yr)];   units="years")

        # State variables
        _defvar(ds, "ocn_temp",     Float32, ("lon","lat","z"),  state.T;        units="degC")
        _defvar(ds, "ocn_salt",     Float32, ("lon","lat","z"),  state.S;        units="psu")
        _defvar(ds, "ocn_u",        Float32, ("lon","lat","z"),  state.u;        units="m/s")
        _defvar(ds, "ocn_v",        Float32, ("lon","lat","z"),  state.v;        units="m/s")
        _defvar(ds, "ocn_w",        Float32, ("lon","lat","zf"), state.w;        units="m/s")
        _defvar(ds, "ocn_psi_bt",   Float32, ("lon","lat"),      state.psi_bt;   units="m2/s")
        _defvar(ds, "ocn_rho",      Float32, ("lon","lat","z"),  state.rho;      units="kg/m3")
        _defvar(ds, "ocn_pressure", Float32, ("lon","lat","z"),  state.pressure; units="bar")
    end
    return filename
end

function _defvar(ds, name, T, dims, data; units="")
    v = NCDatasets.defVar(ds, name, T, dims)
    units != "" && (v.attrib["units"] = units)
    v[:] = data
    return v
end

function _load_ncdatasets()
    try
        return Base.require(Base.PkgId(Base.UUID("85f8d34a-cbdd-5861-8df4-14fed0d494ab"),
                                        "NCDatasets"))
    catch
        error("NCDatasets.jl is required for NetCDF output. Add it with `Pkg.add(\"NCDatasets\")`.")
    end
end

end # module Output
