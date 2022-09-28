# WVZAnalysis

## Words for Devs
1. Main looper logic is in `src/mainlooper.jl`
2. Cuts logic for different channel are in `src/{ZZZ, WZZ, WWZ}.jl`
3. For histogram we use https://github.com/Moelf/FHist.jl
4. I would recomment to use https://github.com/andyferris/Dictionaries.jl for in-memory "hadd" when combining results from different files.
5. all public / user-facing functions are in `analysis_utils.jl`
6. all supporting for printing reports and plotting are in `reporting_utils.jl`

## Install the dependencies:
### Already cloned
If you want the source files to be somewhere else, you can manually `git clone` this package. Then navigate
to that loacation:
```
julia --project=.
] instantiate
```
(the dot `.` just means current directory, as usual)


## Helpful setup cell:
```julia
using Pkg, WVZAnalysis

Pkg.activate(pathof(WVZAnalysis) |> dirname |> dirname)

using ProgressMeter, PrettyTables, UnROOT, FHist, JSON3, CairoMakie,
      ThreadsX, Measurements, Arrow, Serialization
```

## Example Usage:

### Shape syst:
```julia
Hs_up = shapesys("Signal", "EG_SCALE_ALL__1up"; NN_hist=true)
Hs_nominal = sfsys("Signal"; NN_hist=true)
Hs_down = shapesys("Signal", "EG_SCALE_ALL__1down"; NN_hist=true;

Hs = merge(Hs_up, Hs_nominal, Hs_down)
```

### Filling Arrow files
```julia
for tag in WVZAnalysis.ALL_TAGS
    fname = "/data/jiling/WVZ/v2.3_arrow/$tag.arrow"
    data = arrow_making(tag)
    Arrow.write(fname, Dict(pairs(data)))
end
```

### Making Significance table
```julia
# use `recreate=false` if you want to read from `/data/jiling` on af.uchicago
julia> M = significance_table(; recreate=true);

julia> print_sigtable(M)
┌───────────────┬────────────────┬────────────────┬──────────────┐
│               │     SF-inZ     │     SF-noZ     │      DF      │
├───────────────┼────────────────┼────────────────┼──────────────┤
│    Signal     │  11.73 ± 3.43  │  9.31 ± 3.05   │ 10.73 ± 3.28 │
│      ZZ       │ 1775.6 ± 42.14 │ 469.06 ± 21.66 │ 19.78 ± 4.45 │
│     Zjets     │  -0.02 ± 0.14  │   2.6 ± 1.61   │ 6.47 ± 2.54  │
│    Zgamma     │   0.0 ± 0.0    │   0.0 ± 0.0    │  0.3 ± 0.55  │
│     ttbar     │   0.0 ± 0.0    │  0.63 ± 0.79   │ 0.28 ± 0.53  │
│      WZ       │   0.36 ± 0.6   │  1.79 ± 1.34   │  2.24 ± 1.5  │
│      tZ       │  0.01 ± 0.11   │  0.07 ± 0.27   │ 0.06 ± 0.25  │
│      ttZ      │  1.24 ± 1.11   │  4.71 ± 2.17   │  5.74 ± 2.4  │
│      tWZ      │  0.57 ± 0.75   │  2.16 ± 1.47   │  2.5 ± 1.58  │
│      VBS      │  11.71 ± 3.42  │   6.4 ± 2.53   │ 0.18 ± 0.43  │
│    Others     │  0.06 ± 0.24   │   0.4 ± 0.63   │ 0.56 ± 0.75  │
├───────────────┼────────────────┼────────────────┼──────────────┤
│   Bkg Tot.    │ 1789.53 ± 42.3 │ 487.81 ± 22.09 │ 38.11 ± 6.17 │
├───────────────┼────────────────┼────────────────┼──────────────┤
│ Significance  │  0.28 ± 0.08   │  0.42 ± 0.14   │ 1.74 ± 0.55  │
│ Combined Sig. │   NaN ± 0.0    │  1.81 ± 0.53   │  NaN ± 0.0   │
└───────────────┴────────────────┴────────────────┴──────────────┘
```
