# Integration test placeholder for multi-century spinup.
#
# These tests require the reference fixtures (test/fixtures/*.nc) which are
# generated from the CLIMBER-X Fortran model.  They are skipped if fixtures
# are absent (CI without fixture data).
#
# To run against fixtures, place the NetCDF files in test/fixtures/ and run:
#   julia --project test/runtests.jl

using Test

const FIXTURE_DIR = joinpath(@__DIR__, "..", "fixtures")

@testset "Integration / spinup (requires fixtures)" begin

    fixture_10yr  = joinpath(FIXTURE_DIR, "snapshot_year0010.nc")
    fixture_100yr = joinpath(FIXTURE_DIR, "snapshot_year0100.nc")

    if !isfile(fixture_10yr)
        @warn "Fixture $(fixture_10yr) not found; skipping integration tests."
        @test_skip "10-year spinup test skipped (no fixture)"
        return
    end

    @testset "10-year spinup density check" begin
        @test_skip "Fixture-based integration tests not yet implemented; place NC files in test/fixtures/"
    end

    @testset "100-year spinup T/S check" begin
        @test_skip "Fixture-based integration tests not yet implemented"
    end

    @testset "AMOC at 26°N within 3 Sv" begin
        @test_skip "Fixture-based integration tests not yet implemented"
    end
end
