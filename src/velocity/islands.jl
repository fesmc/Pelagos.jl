# Island integral constraints for the barotropic streamfunction solver.
#
# In a multiply-connected ocean domain, the barotropic vorticity equation (an
# elliptic PDE for ψ) has one undetermined integration constant per island.
# These constants are fixed by requiring zero net normal flow around each island
# perimeter:  ∮ (∂ψ/∂n) dl = 0.
#
# This module provides:
#   1. Island detection from the land-sea mask.
#   2. Assembly of the island-integral constraint rows for the sparse system.
#
# Reference: Fortran source ocn_islands.f90 in CLIMBER-X.

module Islands

using SparseArrays

export detect_islands, IslandInfo

"""
    IslandInfo

Holds connectivity data for one island.

Fields:
- `id`           : island index (1-based)
- `perimeter`    : vector of (i,j) tuples marking ocean cells adjacent to this island
- `interior_mask`: Bool matrix marking the island's land interior cells
"""
struct IslandInfo
    id           :: Int
    perimeter    :: Vector{Tuple{Int,Int}}   # ocean cells bordering island
    interior_mask:: Matrix{Bool}             # land interior of this island
end

"""
    detect_islands(ocean_mask) -> Vector{IslandInfo}

Detect islands (connected land regions not touching the domain boundary) from
the binary ocean mask. Returns one `IslandInfo` per island.

The main continent + domain boundary are excluded — only interior land cells
enclosed by ocean on all sides count as islands.

`ocean_mask`: Bool (nlon, nlat), true = ocean cell.
"""
function detect_islands(ocean_mask::Matrix{Bool})::Vector{IslandInfo}
    nlon, nlat = size(ocean_mask)
    land_mask  = .!ocean_mask

    # Label connected land regions using a simple flood-fill
    labels = zeros(Int, nlon, nlat)
    current_label = 0

    for j in 1:nlat, i in 1:nlon
        if land_mask[i, j] && labels[i, j] == 0
            current_label += 1
            _flood_fill!(labels, land_mask, i, j, current_label, nlon, nlat)
        end
    end

    n_regions = current_label
    islands = IslandInfo[]

    for lbl in 1:n_regions
        # Determine if this region touches the domain boundary
        touches_boundary = false
        for j in 1:nlat, i in 1:nlon
            if labels[i, j] == lbl
                if i == 1 || i == nlon || j == 1 || j == nlat
                    touches_boundary = true
                    break
                end
            end
        end
        touches_boundary && continue   # skip main continent

        # Build interior mask and find perimeter ocean cells
        interior = labels .== lbl
        perimeter = Tuple{Int,Int}[]
        for j in 1:nlat, i in 1:nlon
            interior[i, j] || continue
            for (di, dj) in ((-1,0),(1,0),(0,-1),(0,1))
                ni = mod1(i + di, nlon)
                nj = j + dj
                (1 ≤ nj ≤ nlat) || continue
                if ocean_mask[ni, nj]
                    push!(perimeter, (ni, nj))
                end
            end
        end
        unique!(perimeter)
        push!(islands, IslandInfo(length(islands)+1, perimeter, interior))
    end
    return islands
end

function _flood_fill!(labels::Matrix{Int}, land::AbstractMatrix{Bool},
                      i0::Int, j0::Int, lbl::Int,
                      nlon::Int, nlat::Int)
    stack = [(i0, j0)]
    while !isempty(stack)
        i, j = pop!(stack)
        (1 ≤ i ≤ nlon && 1 ≤ j ≤ nlat) || continue
        (!land[i, j] || labels[i, j] != 0) && continue
        labels[i, j] = lbl
        # 4-connected neighbours; longitude is periodic
        push!(stack, (mod1(i-1, nlon), j), (mod1(i+1, nlon), j),
                     (i, j-1),              (i, j+1))
    end
end

end # module Islands
