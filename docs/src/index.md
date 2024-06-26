```@raw html
---
layout: home

hero:
  name: "WVZAnalysis.jl"
  text: Run 2 VVZ analysis (4-lepton channel)
  tagline: "ANA-STDM-2020-08"
  image:
    src: /logo.png
    alt: DocumenterVitepress
  actions:
    - theme: brand
      text: Walkthrough
      link: /walkthrough
    - theme: alt
      text: APIs
      link: /internalapis
    - theme: alt
      text: ATLAS-Glance
      link: "https://atlas-glance.cern.ch/atlas/analysis/analyses/details.php?ref_code=ANA-STDM-2020-08"

---
```

## The Big Picture

- like any Julia package, the source files are under [/src](https://github.com/Moelf/WVZAnalysis.jl/tree/master/src)
- meta data files such as tag -> DSID mapping are under [/config](https://github.com/Moelf/WVZAnalysis.jl/tree/master/config)
- the main looper is called `main_looper()` and all the real actions are in there
- the atomic unit of runnable task is called `AnalysisTask`

## A typical workflow
Say we're trying to produce BDT score histogram for process `ZZ`:

1. call `prep_tasks("ZZ")` (this returns a vector of `AnalysisTasks`)
2. run `main_looper` over each of the tasks, you can use threading, or distributed parallelism.
3. each `main_looper` would would return a dictionary of histograms
4. `merge` all the result (`reduce(mergewith(+), results_from_3)`)
5. profit


## Public-ish APIs
These are the types and functions all a user would need, 
in order to do pretty much anything.

The packages that are good to be familiar with: 

- [UnROOT.jl](https://github.com/JuliaHEP/UnROOT.jl)
- [FHist.jl](https://github.com/Moelf/FHist.jl)


These functions/types and the direct callees are more documented than the rest of the library:

### Analysis and data dumping
```@docs
WVZAnalysis.AnalysisTask
WVZAnalysis.prep_tasks
WVZAnalysis.main_looper(t::AnalysisTask)
WVZAnalysis.arrow_main
```

### Result reporting
```@docs
WVZReportExt.significance_table
WVZReportExt.print_sigtable
```
