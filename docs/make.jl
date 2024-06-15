using Documenter
using DocumenterVitepress

using WVZAnalysis, WVZReportExt

repopath = "github.com/Moelf/WVZAnalysis.jl"

makedocs(;
         modules=[WVZAnalysis, WVZReportExt],
         format=DocumenterVitepress.MarkdownVitepress(; repo = repopath, devbranch = "main", devurl = "dev"),
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
           repo=repopath
          )
