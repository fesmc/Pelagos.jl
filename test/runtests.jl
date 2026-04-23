using Pelagos
using Test

@testset "Pelagos.jl" begin
    @testset "EOS" begin
        include("eos/test_UNESCO.jl")
        include("eos/test_pressure.jl")
    end
    @testset "Velocity" begin
        include("velocity/test_islands.jl")
        include("velocity/test_baroclinic.jl")
        include("velocity/test_barotropic.jl")
        include("velocity/test_continuity.jl")
    end
    @testset "Tracers" begin
        include("tracers/test_diffusion.jl")
        include("tracers/test_forcing.jl")
    end
    @testset "Integration" begin
        include("integration/test_spinup.jl")
    end
end
