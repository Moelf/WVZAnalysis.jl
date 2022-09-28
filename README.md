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
Hs_up = WVZAnalysis.shapesys("Signal", "EG_SCALE_ALL__1up")
Hs_nominal = WVZAnalysis.sfsys("Signal")
Hs_down = WVZAnalysis.shapesys("Signal", "EG_SCALE_ALL__1down");
```

### Filling Arrow files
```julia
for tag in ALL_TAGS
    fname = "/data/jiling/WVZ/v2.3_arrow/$tag.arrow"
    data = WVZAnalysis.arrow_making(tag)
    Arrow.write(fname, Dict(pairs(data)))
end
```

### Making Significance table
```julia
function significance_table()
    proc_names = ("Signal", "ZZ", "Zjets", "Zgamma", "ttbar", "WZ", "tZ", "ttZ", "tWZ", "VBS", "Others")
    M = mapreduce(vcat, proc_names) do tag
        ## re-make
        res = WVZAnalysis.sfsys(tag; NN_hist=true)
        serialize("/data/jiling/WVZ/v2.3_hists/$(tag).jlserialize", res)
        ## load from serialization
        # res = deserialize("/data/jiling/WVZ/v2.3_hists/$(tag).jlserialize")
        N = nbins(res[:DF__NN__NOMINAL])
        hists = rebin.([res[:SFinZ__NN__NOMINAL], res[:SFnoZ__NN__NOMINAL], res[:DF__NN__NOMINAL]], N)
        permutedims(hists) 
    end
    proc_names, M
end

function significance_table(proc_names, M)
    body = @. integral(M) ± only(binerrors(M))
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

### 
