# WVZAnalysis

## Words for Devs
1. Main looper logic is in `src/mainlooper.jl`
2. Cuts logic for different channel are in `src/ZZZ WZZ WWZ`
3. For histogram we use https://github.com/Moelf/FHist.jl
4. I would recomment to use https://github.com/andyferris/Dictionaries.jl for in-memory "hadd" when combining results from different files.

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

using ProgressMeter, PrettyTables, UnROOT, FHist, JSON3, CairoMakie, ThreadsX, Arrow, Measurements

BKG_TAGS = ("ZZ", "Zjets", "Zgamma", "ttbar", "WZ", "tZ", "ttZ", "tWZ", "VBS", "VH", "Others")
ALL_TAGS = ("Signal", BKG_TAGS...)
```

## Example Usage:

### Shape syst:
```julia
Hs_up = WVZAnalysis.shapesys("Signal", "tree_EG_SCALE_ALL__1up")
Hs_nominal = WVZAnalysis.sfsys("Signal")
Hs_down = WVZAnalysis.shapesys("Signal", "tree_EG_SCALE_ALL__1down");
```

### Filling Arrow files
```julia
for tag in ALL_TAGS
    fname = "/data/jiling/WVZ/v2.2_arrow/$tag.arrow"
    data = WVZAnalysis.arrow_making(tag)
    Arrow.write(fname, Dict(pairs(data)))
end
```

### Making Significance table
```julia
function significance_table()
    proc_names = ALL_TAGS
    M = mapreduce(vcat, proc_names) do tag
        res = WVZAnalysis.sfsys(tag)
        hists = [res[:NN_inZ__NOMINAL], res[:NN_noZ__NOMINAL], res[:NN_DF__NOMINAL]]
        permutedims(hists) # make each sample a row
    end
    proc_names, M
end

function significance_table(proc_names, M)
    body = @. integral(M) ± (only(binerrors(M)))
    total_sig = body[1:1, :] #first row
    total_bkg = sum(body[2:end, :]; dims = 1) #2:end row
    naive_significance = @. total_sig / sqrt(total_bkg)
    combined_sig = sqrt(sum(abs2, naive_significance))

    full_nums = [
        body
        total_bkg
        naive_significance
        [NaN combined_sig NaN]
    ]

    full_body = [collect(proc_names)
                 "Bkg Tot."
                 "Significance"
                 "Combined Sig.";; full_nums
                 ];
end

proc_names, M = significance_table();

full_table = significance_table(proc_names, M);

fmt = (v, i, j) -> v isa Number ? "$(round(Measurements.value(v); digits=2)) ± $(round(Measurements.uncertainty(v); digits=2))" : v
pretty_table(
    full_table;
    header = ["", "SF-inZ", "SF-noZ", "DF"], 
    formatters = fmt,
    body_hlines = [size(full_table, 1) - 3, size(full_table, 1) - 2],
    highlighters = hl_row([1], crayon"bold"), 
    alignment=:c,
    #backend = Val(:html)
)
┌───────────────┬────────────────┬───────────────┬──────────────┐
│               │     SF-inZ     │    SF-noZ     │      DF      │
├───────────────┼────────────────┼───────────────┼──────────────┤
│    Signal     │  13.35 ± 0.07  │ 10.69 ± 0.11  │ 12.26 ± 0.15 │
│      ZZ       │ 2049.2 ± 6.13  │ 541.18 ± 2.59 │ 23.16 ± 0.48 │
│     Zjets     │  0.64 ± 0.34   │  3.67 ± 2.22  │ 6.72 ± 5.16  │
│    Zgamma     │   0.0 ± 0.0    │  0.31 ± 0.31  │  0.2 ± 0.31  │
│     ttbar     │  0.00 ± 0.00   │  0.96 ± 0.21  │ 0.55 ± 0.14  │
│      WZ       │  0.68 ± 0.14   │  3.28 ± 0.31  │ 3.72 ± 0.35  │
│      tZ       │  0.04 ± 0.02   │  0.17 ± 0.05  │ 0.11 ± 0.03  │
│      ttZ      │  1.46 ± 0.09   │  5.41 ± 0.17  │ 6.61 ± 0.19  │
│      tWZ      │  0.71 ± 0.57   │  2.39 ± 0.24  │ 3.09 ± 0.26  │
│      VBS      │  13.22 ± 0.1   │  7.12 ± 0.08  │ 0.23 ± 0.01  │
│      VH       │  1.29 ± 0.71   │  5.9 ± 1.41   │ 6.32 ± 1.35  │
│    Others     │  0.07 ± 0.01   │  0.6 ± 0.2    │  0.6 ± 0.08  │
├───────────────┼────────────────┼───────────────┼──────────────┤
│   Bkg Tot.    │ 2067.37 ± 7.88 │ 570.99 ± 7.79 │ 51.33 ± 8.37 │
├───────────────┼────────────────┼───────────────┼──────────────┤
│ Significance  │   0.29 ± 0.0   │  0.45 ± 0.0   │ 1.71 ± 0.07  │
│ Combined Sig. │   NaN ± 0.0    │  1.6 ± 0.06   │  NaN ± 0.0   │
└───────────────┴────────────────┴───────────────┴──────────────┘
```
