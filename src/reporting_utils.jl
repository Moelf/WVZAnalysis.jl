using PrettyTables, Serialization, Measurements

"""
    rebinscan(S, B; atleast=1, from=:right, by = (s,b) -> s/sqrt(b))

Take two `Hist1D`, try to find the best bin-edges for rebinning, given these two conditions:

- a metric to maximize, `s` and `b` are signal counts and background counts in each new bin
- `atleast` specifcy the minimal number of `s+b` in each newly formed bin

`from` keyword argument can be used to specifcy the direction of the scan, defautls to `:right` assuming
the high-er end of the histogram is signal-like and one of the histograms is monotonic.

The function returns the values of new bin-edges including both ends (of course, the ends are identical to before).

The metric can also be replaced by the more accurate `sqrt(2*((s + b) * log(1 + s/b) - s ))`.
"""
function rebinscan(S, B; atleast=1, from=:right, by = (s,b) -> s/sqrt(b))
    _binedges = binedges(S)
    _binedges == binedges(B) || error("Bin edges aren't compatible")
    Scounts = bincounts(S)
    Bcounts = bincounts(B)
    if from == :right
        Scounts = reverse(Scounts)
        BCounts = reverse(Bcounts)
    end

    tempS = tempB = 0.0
    score = -Inf
    newEdges = [1] # first edge at the beginning
    for i = 2:lastindex(_binedges) # one more edge than bins
        s = Scounts[i-1]
        b = Bcounts[i-1]
        tempS += s
        tempB += b

        tempS + tempB < atleast && continue # not enough statistics
        newScore = by(tempS, tempB)
        if newScore < score # score starts to drop
            push!(newEdges, i-1)
            tempS = s
            tempB = b
            score = -Inf
        end
        score = newScore
    end
    push!(newEdges, lastindex(_binedges))
    return _binedges[newEdges]
end

function _significance(signal, bkg)
    S = integral(signal)
    B = integral(bkg)
    Sig = sqrt(2*((S + B) * log(1 + S/B) - S))
    dS = sqrt(sum(signal.sumw2))
    dB = sqrt(sum(bkg.sumw2))
    dSigdS = log(1 + S/B) / Sig
    dSigdB = (log(1 + S/B) - S/B) / Sig
    err = sqrt((dSigdS * dS)^2 + (dSigdB * dB)^2)
    return Sig ± err
end

function significance_matrix(; recreate, mapper=map)
    Ms = map(ALL_TAGS) do tag
        deserialize(joinpath(ANALYSIS_DIR[],"$(tag).jlserialize"))
    end
    significance_matrix(Ms)
end

function significance_matrix(Ms)
    res = mapreduce(vcat, Ms) do M
        N = nbins(M[:DF__BDT__NOMINAL])
        hists = rebin.([M[:SFinZ__BDT__NOMINAL], M[:SFnoZ__BDT__NOMINAL], M[:DF__BDT__NOMINAL]], N)
        N = nbins(M[:ZZCR__Njet__NOMINAL])
        push!(hists, rebin(M[:ZZCR__Njet__NOMINAL], N))
        N = nbins(M[:ttZCR__Njet__NOMINAL])
        push!(hists, rebin(M[:ttZCR__Njet__NOMINAL], N))

        permutedims(hists) 
    end
    res
end

"""
    significance_table(; recreate=false)
    significance_table(significance_matrix(); recreate=false)

# Examples

```julia
julia> M = significance_table()
14×6 Matrix{Any}:
 "Signal"         10.66±0.067   …    1.072±0.017     2.151±0.037
 "ZZ"            1219.9±5.4          555.7±2.4       69.19±0.74
 "Zjets"         -0.019±0.13       -0.0004±0.0004     0.63±0.39
 "Zgamma"           0.0±0.0            0.0±0.0         0.0±0.0
 "ttbar"            0.0±0.0            0.0±0.0        0.43±0.13
 "WZ"             0.358±0.096   …   0.0044±0.0032    0.295±0.061
 "tZ"             0.012±0.012          0.0±0.0       0.179±0.043
 "ttZ"            1.231±0.08        0.0111±0.0074    73.62±0.61
 "tWZ"             0.56±0.11        0.0095±0.0095    14.71±0.56
 "VBS"           10.231±0.085        1.478±0.033     0.833±0.024
 "VH"              1.29±0.71    …      0.0±0.0         0.0±0.0
 "Others"        0.0568±0.0088      0.0021±0.0017     5.02±0.088
 "Bkg Tot."      1233.6±5.4          557.2±2.4       164.9±1.2
 "Significance"  0.3031±0.002       0.0454±0.00072  0.1672±0.0029
```
"""
function significance_table(; recreate=false)
    body = significance_matrix(; recreate)
    significance_table(body)
end

function significance_table(body::Matrix; recreate=false)
    total_sig = body[1:1, :] #first row
    total_bkg = mapreduce(sum, hcat, eachcol(body[2:end, :])) #2:end row
    sig_errors = _significance.(total_sig, total_bkg)
    combined_sig = sqrt(sum(x[1]^2 for x in sig_errors))

    full_nums = [
        @. integral(body) ± only(binerrors(body))
        @. integral(total_bkg) ± only(binerrors(total_bkg))
        sig_errors
    ]

    full_body = [collect(ALL_TAGS)
                 "Bkg Tot."
                 "Significance" ;; full_nums
                 ];
end

const sigtable_fmt = (v, i, j) -> v isa Number ? "$(round(Measurements.value(v); digits=2)) ± $(round(Measurements.uncertainty(v); digits=2))" : v

"""
    print_sigtable(full_table; io=stdout)

Takes the output of [`significance_table`](@ref) and pretty print it:

# Example

```julia
julia> M = significance_table(; recreate=true);

julia> print_sigtable(M)
┌──────────────┬────────────────┬───────────────┬──────────────┬──────────────┬───────────────┐
│              │     SF-inZ     │    SF-noZ     │      DF      │    CR-ZZ     │    CR-ttZ     │
├──────────────┼────────────────┼───────────────┼──────────────┼──────────────┼───────────────┤
│    Signal    │  10.66 ± 0.07  │  9.31 ± 0.1   │ 10.73 ± 0.14 │ 1.07 ± 0.02  │  2.15 ± 0.04  │
│      ZZ      │ 1219.92 ± 5.38 │ 469.06 ± 2.44 │ 19.78 ± 0.45 │ 555.68 ± 2.4 │ 69.19 ± 0.74  │
│    Zjets     │  -0.02 ± 0.13  │  2.6 ± 2.22   │ 6.47 ± 5.51  │  -0.0 ± 0.0  │  0.63 ± 0.39  │
│    Zgamma    │   0.0 ± 0.0    │   0.0 ± 0.0   │  0.3 ± 0.29  │  0.0 ± 0.0   │   0.0 ± 0.0   │
│    ttbar     │   0.0 ± 0.0    │  0.63 ± 0.18  │  0.28 ± 0.1  │  0.0 ± 0.0   │  0.43 ± 0.13  │
│      WZ      │   0.36 ± 0.1   │  1.79 ± 0.23  │ 2.24 ± 0.29  │  0.0 ± 0.0   │  0.29 ± 0.06  │
│      tZ      │  0.01 ± 0.01   │  0.07 ± 0.03  │ 0.06 ± 0.02  │  0.0 ± 0.0   │  0.18 ± 0.04  │
│     ttZ      │  1.23 ± 0.08   │  4.71 ± 0.16  │ 5.74 ± 0.18  │ 0.01 ± 0.01  │ 73.62 ± 0.61  │
│     tWZ      │  0.56 ± 0.11   │  2.16 ± 0.23  │  2.5 ± 0.24  │ 0.01 ± 0.01  │ 14.71 ± 0.56  │
│     VBS      │  10.23 ± 0.09  │  6.4 ± 0.08   │ 0.18 ± 0.01  │ 1.48 ± 0.03  │  0.83 ± 0.02  │
│      VH      │  1.29 ± 0.71   │  5.76 ± 1.4   │ 5.77 ± 1.29  │  0.0 ± 0.0   │   0.0 ± 0.0   │
├──────────────┼────────────────┼───────────────┼──────────────┼──────────────┼───────────────┤
│    Others    │  0.06 ± 0.01   │  0.4 ± 0.13   │ 0.56 ± 0.08  │  0.0 ± 0.0   │  5.02 ± 0.09  │
├──────────────┼────────────────┼───────────────┼──────────────┼──────────────┼───────────────┤
│   Bkg Tot.   │ 1233.64 ± 5.43 │ 493.58 ± 3.61 │ 43.89 ± 5.7  │ 557.19 ± 2.4 │ 164.91 ± 1.19 │
│ Significance │   0.3 ± 0.0    │  0.42 ± 0.0   │  1.56 ± 0.1  │  0.05 ± 0.0  │  0.17 ± 0.0   │
└──────────────┴────────────────┴───────────────┴──────────────┴──────────────┴───────────────┘
```
"""
print_sigtable(full_table; io=stdout) = pretty_table(io, 
    full_table;
    header = ["", "SF-inZ", "SF-noZ", "DF", "CR-ZZ", "CR-ttZ"],
    formatters = sigtable_fmt,
    body_hlines = [size(full_table, 1) - 3, size(full_table, 1) - 2],
    highlighters = hl_row([1], crayon"bold"), 
    alignment=:c,
)
