# Checkpoint read/write for Pelagos.jl.
#
# Restarts store the full ocean state (T, S) and the Oceananigans model state
# to a NetCDF file.  The velocity fields (u, v, w, ψ) are diagnosed each
# timestep and do not need to be checkpointed.

module Restarts

export write_restart, read_restart

"""
    write_restart(filename, T, S, time_yr)

Write a restart file containing tracer fields and the current model time.

# Arguments
- `filename`: path to output .nc file
- `T`       : temperature (nlon, nlat, nz), °C
- `S`       : salinity (nlon, nlat, nz), psu
- `time_yr` : current model time, years
"""
function write_restart(filename::String,
                       T       ::AbstractArray{Float64,3},
                       S       ::AbstractArray{Float64,3},
                       time_yr ::Float64)
    NCDatasets = _load_ncdatasets()
    nlon, nlat, nz = size(T)
    NCDatasets.Dataset(filename, "c") do ds
        NCDatasets.defDim(ds, "nlon", nlon)
        NCDatasets.defDim(ds, "nlat", nlat)
        NCDatasets.defDim(ds, "nz",   nz)

        v_time = NCDatasets.defVar(ds, "time", Float64, ())
        v_time.attrib["units"] = "years"
        v_time[] = time_yr

        v_T = NCDatasets.defVar(ds, "T", Float64, ("nlon","nlat","nz"))
        v_T.attrib["units"] = "degC"
        v_T[:] = T

        v_S = NCDatasets.defVar(ds, "S", Float64, ("nlon","nlat","nz"))
        v_S.attrib["units"] = "psu"
        v_S[:] = S
    end
    return filename
end

"""
    read_restart(filename) -> (T, S, time_yr)

Read a restart file.  Returns tracer arrays and model time.
"""
function read_restart(filename::String)
    NCDatasets = _load_ncdatasets()
    NCDatasets.Dataset(filename, "r") do ds
        T       = Float64.(ds["T"][:])
        S       = Float64.(ds["S"][:])
        time_yr = Float64(ds["time"][])
        return (T, S, time_yr)
    end
end

function _load_ncdatasets()
    try
        return Base.require(Base.PkgId(Base.UUID("85f8d34a-cbdd-5861-8df4-14fed0d494ab"),
                                        "NCDatasets"))
    catch
        error("NCDatasets.jl is required for restart I/O. Add it with `Pkg.add(\"NCDatasets\")`.")
    end
end

end # module Restarts
