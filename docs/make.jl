using Documenter
using Pelagos

makedocs(;
    modules   = [Pelagos],
    sitename  = "Pelagos.jl",
    format    = Documenter.HTML(;
        prettyurls    = get(ENV, "CI", "false") == "true",
        canonical     = "https://fesmc.github.io/Pelagos.jl/",
        edit_link     = "main",
        assets        = String[],
    ),
    pages = [
        "Home"        => "index.md",
        "Physics"     => "physics.md",
        "Architecture"=> "architecture.md",
        "API"         => [
            "Equation of State" => "api/eos.md",
            "Velocity Solver"   => "api/velocity.md",
            "Tracer Model"      => "api/tracers.md",
            "Forcing"           => "api/forcing.md",
            "Grid"              => "api/grid.md",
            "I/O"               => "api/io.md",
        ],
        "Development" => "development.md",
    ],
    checkdocs = :none,
)

deploydocs(;
    repo      = "github.com/fesmc/Pelagos.jl.git",
    devbranch = "main",
)
