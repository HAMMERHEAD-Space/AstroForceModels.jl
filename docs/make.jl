using Documenter
using AstroForceModels

makedocs(;
    modules=[AstroForceModels],
    format=Documenter.HTML(;
        prettyurls=(!("local" in ARGS)), highlights=["yaml"], ansicolor=true
    ),
    sitename="AstroForceModels.jl",
    authors="Jordan Murphy",
    pages=[
        "Home" => "index.md",
        "Usage" => "man/usage.md",
        "Force Models" => Any[
            "force_models/gravity.md",
            "force_models/drag.md",
            "force_models/solar_radiation_pressure.md",
            "force_models/third_body.md",
            "force_models/relativity.md",
        ],
        "API Reference" => "man/api.md",
        "Library" => "lib/library.md",
    ],
)

deploydocs(; repo="github.com/jmurphy6895/AstroForceModels.jl.git", target="build")
