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

const BKG_TAGS = ("ZZ", "Zjets", "Zgamma", "ttbar", "WZ", "tZ", "ttZ", "tWZ", "VBS", "VH", "Others")
const ALL_TAGS = ("Signal", BKG_TAGS)
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
│    Signal     │  13.68 ± 0.08  │ 11.01 ± 0.12  │ 12.42 ± 0.15 │
│      ZZ       │ 2110.05 ± 6.33 │ 604.98 ± 2.82 │ 32.4 ± 0.57  │
│     Zjets     │  0.64 ± 0.34   │  2.7 ± 1.11   │ 14.05 ± 5.62 │
│    Zgamma     │   0.0 ± 0.0    │  0.38 ± 0.33  │  0.33 ± 0.3  │
│     ttbar     │  0.07 ± 0.05   │  0.84 ± 0.2   │ 0.76 ± 0.17  │
│      WZ       │  0.52 ± 0.13   │  3.15 ± 0.31  │ 3.38 ± 0.34  │
│      tZ       │  0.02 ± 0.01   │  0.15 ± 0.06  │ 0.12 ± 0.04  │
│      ttZ      │  1.59 ± 0.09   │  6.03 ± 0.18  │  7.29 ± 0.2  │
│      tWZ      │  0.78 ± 0.13   │  2.68 ± 0.26  │  3.2 ± 0.27  │
│      VBS      │   13.6 ± 0.1   │  7.95 ± 0.08  │ 0.29 ± 0.01  │
│      VH       │  2.11 ± 0.91   │  6.75 ± 1.51  │ 4.83 ± 1.14  │
│    Others     │  0.07 ± 0.01   │  0.5 ± 0.14   │ 0.58 ± 0.08  │
├───────────────┼────────────────┼───────────────┼──────────────┤
│   Bkg Tot.    │ 2129.44 ± 6.41 │ 636.11 ± 3.45 │ 67.22 ± 5.79 │
├───────────────┼────────────────┼───────────────┼──────────────┤
│ Significance  │   0.3 ± 0.0    │  0.44 ± 0.0   │ 1.52 ± 0.07  │
│ Combined Sig. │   NaN ± 0.0    │  1.6 ± 0.06   │  NaN ± 0.0   │
└───────────────┴────────────────┴───────────────┴──────────────┘
```
