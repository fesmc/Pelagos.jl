# CLIMBER-X compatible Oceananigans grid construction.
#
# The ocean grid is 72×36×23 at 5°×5° resolution.  Layer interfaces are read
# directly from the `zw` variable in an ocn_restart.nc file so that the
# Oceananigans grid matches the CLIMBER-X vertical levels exactly.
#
# Bathymetry is derived from `k1_pot` (bottom layer index) and `f_ocn`
# (ocean fraction) in the same restart file, and applied via
# Oceananigans' GridFittedBottom immersed boundary.

module GridSetup

using Oceananigans
using Oceananigans.Grids: LatitudeLongitudeGrid
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom
using NCDatasets

export build_climberx_grid

const CLIMBERX_NX = 72
const CLIMBERX_NY = 36
const CLIMBERX_NZ = 23

"""
    build_climberx_grid(arch=CPU(); restart_file) -> ImmersedBoundaryGrid

Build an Oceananigans ImmersedBoundaryGrid matching the CLIMBER-X ocean model:
- 72 × 36 × 23, 5° × 5° resolution
- Vertical interfaces read from `zw` in the restart NetCDF (ordered bottom→top)
- Bathymetry derived from `k1_pot` (bottom layer index) and `f_ocn` (ocean mask)

Land cells have `bottom_height = 0` (sealing the full water column from above).
Ocean cells have `bottom_height = zw[k1_pot[i,j]]` (negative depth in metres).
"""
function build_climberx_grid(arch = CPU(); restart_file::String)
    isfile(restart_file) || error("restart file not found: $restart_file")

    zw, bottom_height = NCDatasets.Dataset(restart_file, "r") do ds
        zw    = Float64.(ds["zw"][:])              # (24,) interfaces, m, bottom→top (negative)
        k1    = reshape(Int.(ds["k1_pot"][:]), CLIMBERX_NX, CLIMBERX_NY)
        f_ocn = reshape(Float64.(ds["f_ocn"][:]), CLIMBERX_NX, CLIMBERX_NY)

        bottom_height = zeros(CLIMBERX_NX, CLIMBERX_NY)
        @inbounds for j in 1:CLIMBERX_NY, i in 1:CLIMBERX_NX
            if f_ocn[i, j] > 0.5
                # zw[k1[i,j]] is the interface below the deepest active layer
                bottom_height[i, j] = zw[k1[i, j]]   # negative (below sea surface)
            end
            # land cells: bottom_height stays 0.0 (seals the column from above)
        end
        zw, bottom_height
    end

    underlying_grid = LatitudeLongitudeGrid(arch;
        size      = (CLIMBERX_NX, CLIMBERX_NY, CLIMBERX_NZ),
        longitude = (-180, 180),
        latitude  = (-90,   90),
        z         = zw,          # 24 interfaces, already bottom→top
        topology  = (Periodic, Bounded, Bounded),
    )

    return ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height))
end

end # module GridSetup
