# Land mask and island detection utilities for Pelagos.jl.
#
# The ocean mask is derived from the (smoothed) bathymetry: a cell is ocean if
# H > 0.  Island detection then identifies connected land regions not touching
# the domain boundary — these require special treatment in the barotropic solver.
#
# Separate masks are maintained for:
#   - Barotropic flow (Bering Strait and Davis Strait: closed)
#   - Tracer transport (Bering Strait and Davis Strait: open)
#
# Island detection is delegated to src/velocity/islands.jl.

module Masks

using ..Islands: IslandInfo, detect_islands

export build_ocean_mask, build_velocity_mask, build_tracer_mask

"""
    build_ocean_mask(H) -> Matrix{Bool}

Derive the ocean mask from the bathymetry depth field.  Cells with H > 0 are ocean.
"""
function build_ocean_mask(H::AbstractMatrix{Float64})::Matrix{Bool}
    return H .> 0.0
end

"""
    build_velocity_mask(ocean_mask; closed_straits=Vector{Tuple{Int,Int}}[]) -> Matrix{Bool}

Build the mask used for barotropic flow.  Optionally close specified grid cells
(e.g., Bering Strait, Davis Strait) by setting them to land.

`closed_straits` is a vector of (i,j) index tuples to close.
"""
function build_velocity_mask(ocean_mask    ::Matrix{Bool};
                             closed_straits::Vector{Tuple{Int,Int}} = Tuple{Int,Int}[])::Matrix{Bool}
    vmask = copy(ocean_mask)
    for (i, j) in closed_straits
        vmask[i, j] = false
    end
    return vmask
end

"""
    build_tracer_mask(ocean_mask) -> Matrix{Bool}

Build the mask used for tracer transport.  Bering and Davis Straits remain open
for tracer exchange even if closed for barotropic flow.  Currently identical to
`ocean_mask`; the distinction is maintained for architectural clarity.
"""
function build_tracer_mask(ocean_mask::Matrix{Bool})::Matrix{Bool}
    return copy(ocean_mask)
end

end # module Masks
