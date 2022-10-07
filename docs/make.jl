using Documenter, WVZAnalysis

makedocs(;
    modules=[WVZAnalysis],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
    ],
    repo="https://github.com/Moelf/WVZAnalysis.jl/blob/{commit}{path}#L{line}",
    sitename="WVZAnalysis.jl",
    authors="Harvard ATLAS",
)

deploydocs(;
    repo="https://github.com/Moelf/WVZAnalysis.jl",
)
