# Tests for island detection from the land-sea mask.

using Test
using Pelagos.Islands: detect_islands, IslandInfo

@testset "Island detection" begin

    @testset "No islands — simple rectangular ocean" begin
        # Land only on the boundary → no interior islands
        nlon, nlat = 8, 6
        mask = fill(true, nlon, nlat)
        mask[:, 1] .= false; mask[:, end] .= false
        mask[1, :] .= false; mask[end, :] .= false

        islands = detect_islands(mask)
        @test isempty(islands)
    end

    @testset "Single interior island" begin
        nlon, nlat = 10, 10
        mask = fill(true, nlon, nlat)
        # Island: a 2×2 block in the interior
        mask[5:6, 5:6] .= false

        islands = detect_islands(mask)
        @test length(islands) == 1
        isl = islands[1]
        @test isl.id == 1
        @test !isempty(isl.perimeter)
        # Perimeter cells must be ocean
        for (i, j) in isl.perimeter
            @test mask[i, j]
        end
    end

    @testset "Two separate islands" begin
        nlon, nlat = 14, 10
        mask = fill(true, nlon, nlat)
        mask[3:4, 4:5]   .= false  # island 1
        mask[10:11, 6:7] .= false  # island 2

        islands = detect_islands(mask)
        @test length(islands) == 2
    end

    @testset "Boundary land is not an island" begin
        nlon, nlat = 8, 8
        mask = fill(true, nlon, nlat)
        # Patch touching the j=1 boundary → should be continent, not island
        mask[4:5, 1:2] .= false

        islands = detect_islands(mask)
        @test isempty(islands)
    end
end
