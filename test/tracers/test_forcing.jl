# Tests for the surface forcing module (virtual salinity flux, global correction).

using Test
using Pelagos.Surface: virtual_salinity_flux!, apply_global_salt_correction!

@testset "Surface forcing" begin

    @testset "Virtual salinity flux" begin
        nlon, nlat = 4, 4
        FW      = fill(1e-8, nlon, nlat)   # 1 cm yr⁻¹ net evap
        S_surf  = fill(35.0, nlon, nlat)
        H_top   = 50.0
        mask    = fill(true, nlon, nlat)
        mask[1,1] = false  # one land cell

        FS = zeros(nlon, nlat)
        virtual_salinity_flux!(FS, FW, S_surf, H_top, mask)

        expected = FW[2,2] * S_surf[2,2] / H_top
        @test isapprox(FS[2, 2], expected, rtol=1e-10)
        # Land cell must be zero
        @test FS[1, 1] == 0.0
    end

    @testset "Global correction gives zero integrated flux" begin
        nlon, nlat = 6, 6
        FS   = randn(nlon, nlat) .* 1e-7
        mask = fill(true, nlon, nlat)
        area = fill(1.0, nlon, nlat)   # equal areas for simplicity

        apply_global_salt_correction!(FS, mask, area)

        total = sum(FS[mask] .* area[mask])
        @test abs(total) < 1e-20
    end
end
