using Documenter, WVZAnalysis

makedocs(;
    modules=[WVZAnalysis],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
        "APIs" => "internalapis.md",
    ],
    repo="https://gitlab.cern.ch/jiling/WVZAnalysis.jl/blob/{commit}{path}#L{line}",
    sitename="WVZAnalysis.jl",
    authors="Jerry Ling and contributors",
)

deploydocs(;
    repo="gitlab.cern.ch/jiling/WVZAnalysis.jl",
)
