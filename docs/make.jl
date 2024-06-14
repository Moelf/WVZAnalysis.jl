using Documenter

using WVZAnalysis, WVZReportExt

makedocs(;
         modules=[WVZAnalysis, WVZReportExt],
         format = Documenter.HTML(
                                  prettyurls = get(ENV, "CI", nothing) == "true",
                                  assets=String[],
                                 ),
         pages=[
                "Introduction" => "index.md",
                "Internal APIs" => "internalapis.md",
                "Step-by-step Walkthrough" => "walkthrough.md",
               ],
         repo="https://github.com/Moelf/WVZAnalysis.jl/blob/{commit}{path}#L{line}",
         sitename="WVZAnalysis.jl",
         authors="Harvard ATLAS",
        )

deploydocs(;
           repo="github.com/Moelf/WVZAnalysis.jl",
          )
