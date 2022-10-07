## The Big Picture

- like any Julia package, the source files are under [/src](./src)
- meta data files such as tag -> DSID mapping are under [/config](./config)
- the main looper is called `main_looper()` and all the real actions are in there
- the atomic unit of runnable task is called `AnalysisTask`

## A typical workflow
Say we're trying to produce BDT score histogram for process `ZZ`:

1. call `prep_tasks("ZZ")` (this returns a vector of `AnalysisTasks`)
2. run `main_looper` over each of the tasks, you can use threading, or distributed
parallelism.
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
AnalysisTask
prep_tasks
main_looper
arrow_making
hist_root
hist_root_pmap
```

### Result reporting
```@docs
significance_table
print_sigtable
```
