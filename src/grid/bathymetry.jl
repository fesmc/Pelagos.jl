# Bathymetry loading and smoothing for Pelagos.jl.
#
# The raw bathymetry is smoothed by a 4-point nearest-neighbour average before
# use (Section 2.2 of Willeit et al. 2022).  This is applied once during grid
# initialisation, not at each timestep.
#
# A topographic slope limiter is also applied to the bathymetry field entering
# the Coriolis term in the barotropic equation (see ocn_grid.f90 in CLIMBER-X).

module Bathymetry

export smooth_bathymetry, load_bathymetry_nc

"""
    smooth_bathymetry(H_raw; n_passes=1) -> Matrix{Float64}

Smooth bathymetry by a 4-point nearest-neighbour kernel.  Longitude is treated
as periodic; northern and southern boundaries use the boundary value as the
neighbour (zero-gradient).

The kernel averages the cell itself and its 4 cardinal neighbours (east, west,
north, south), all with equal weight 1/5:
    H_smooth[i,j] = (H[i,j] + H[E] + H[W] + H[N] + H[S]) / 5

`n_passes` determines how many times the smoothing is applied.  The paper uses
one pass.
"""
function smooth_bathymetry(H_raw::AbstractMatrix{Float64}; n_passes::Int = 1)::Matrix{Float64}
    H = copy(H_raw)
    nlon, nlat = size(H)
    H_tmp = similar(H)

    for _ in 1:n_passes
        @inbounds for j in 1:nlat, i in 1:nlon
            ip1 = mod1(i+1, nlon); im1 = mod1(i-1, nlon)
            jp1 = min(j+1, nlat);  jm1 = max(j-1, 1)
            H_tmp[i,j] = (H[i,j] + H[ip1,j] + H[im1,j] + H[i,jp1] + H[i,jm1]) / 5.0
        end
        H .= H_tmp
    end
    return H
end

"""
    load_bathymetry_nc(filename; lon_var="lon", lat_var="lat", depth_var="depth")
        -> (H, lon, lat)

Load bathymetry from a NetCDF file.  Returns:
- `H`   : depth (nlon, nlat), m, positive downward; land = 0 or NaN
- `lon` : longitude vector, degrees
- `lat` : latitude vector, degrees
"""
function load_bathymetry_nc(filename  ::String;
                            lon_var   ::String = "lon",
                            lat_var   ::String = "lat",
                            depth_var ::String = "depth")
    # NCDatasets loaded at call site to avoid hard dependency at module load
    NCDatasets = Base.require(Base.PkgId(Base.UUID("85f8d34a-cbdd-5861-8df4-14fed0d494ab"),
                                          "NCDatasets"))
    NCDatasets.Dataset(filename) do ds
        lon = Float64.(ds[lon_var][:])
        lat = Float64.(ds[lat_var][:])
        H   = Float64.(ds[depth_var][:,:])
        # Replace missing/NaN with 0 (land)
        H[ismissing.(H) .| isnan.(H)] .= 0.0
        return (H, lon, lat)
    end
end

end # module Bathymetry
