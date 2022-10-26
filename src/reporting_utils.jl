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

function significance(signal, bkg)
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

function significance_matrix(; recreate)
    Ms = map(ALL_TAGS) do tag
        ## re-make
        res = if recreate
            @info tag
            tasks = prep_tasks(tag; NN_hist=true)
            s = ThreadsX.map(main_looper, tasks)
            res = reduce(mergewith(+), s)
            serialize(joinpath(ANALYSIS_DIR[],"$(tag).jlserialize"), res)
            res
        else
            #load from serialization
            deserialize(joinpath(ANALYSIS_DIR[],"$(tag).jlserialize"))
        end
        return res
    end
    significance_matrix(Ms)
end

function significance_matrix(Ms)
    res = mapreduce(vcat, Ms) do M
        N = nbins(M[:DF__NN__NOMINAL])
        hists = rebin.([M[:SFinZ__NN__NOMINAL], M[:SFnoZ__NN__NOMINAL], M[:DF__NN__NOMINAL]], N)
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
15×4 Matrix{Any}:
 "Signal"         11.509±0.068     9.18±0.1     10.14±0.14
 "ZZ"             1709.1±5.5      389.2±2.3      19.1±0.44
 "Zjets"           -0.02±0.13       1.9±1.1       5.9±5.5
 "Zgamma"            0.0±0.0        0.0±0.0       0.3±0.29
 "ttbar"             0.0±0.0       0.63±0.18     0.28±0.1
 "WZ"              0.362±0.096     1.72±0.22     2.13±0.28
 "tZ"                0.0±0.0      0.059±0.025    0.06±0.025
 "ttZ"             1.227±0.08      4.66±0.16     5.68±0.17
 "tWZ"              0.58±0.11      2.15±0.23     2.51±0.24
 "VBS"            11.682±0.092    5.239±0.073   0.182±0.011
 "VH"               1.29±0.71       5.6±1.4       5.4±1.2
 "Others"         0.0584±0.0089    0.38±0.13    0.549±0.082
 "Bkg Tot."       1724.3±5.6      411.5±2.9      42.1±5.7
 "Significance"   0.2768±0.0017  0.4508±0.0053  1.505±0.096
 "Combined Sig."     NaN±0.0      1.595±0.091     NaN±0.0

```
"""
function significance_table(; recreate=false)
    body = significance_matrix(; recreate)
    significance_table(body)
end

function significance_table(body::Matrix; recreate=false)
    total_sig = body[1:1, :] #first row
    total_bkg = mapreduce(sum, hcat, eachcol(body[2:end, :])) #2:end row
    sig_errors = significance.(total_sig, total_bkg)
    combined_sig = sqrt(sum(x[1]^2 for x in sig_errors))

    full_nums = [
        @. integral(body) ± only(binerrors(body))
        @. integral(total_bkg) ± only(binerrors(total_bkg))
        sig_errors
        [NaN combined_sig NaN]
    ]

    full_body = [collect(ALL_TAGS)
                 "Bkg Tot."
                 "Significance"
                 "Combined Sig.";; full_nums
                 ];
end

const sigtable_fmt = (v, i, j) -> v isa Number ? "$(round(Measurements.value(v); digits=2)) ± $(round(Measurements.uncertainty(v); digits=2))" : v

"""
    print_sigtable(full_table)

Takes the output of [`significance_table`](@ref) and pretty print it:

# Example

```julia
julia> M = significance_table(; recreate=true);

julia> print_sigtable(M)

\u250c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u252c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u252c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u252c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2510
\u2502               \u2502     SF-inZ     \u2502    SF-noZ     \u2502      DF      \u2502
\u251c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2524
\u2502    Signal     \u2502  11.73 \u00b1 0.07  \u2502  9.31 \u00b1 0.1   \u2502 10.73 \u00b1 0.14 \u2502
\u2502      ZZ       \u2502 1775.6 \u00b1 5.89  \u2502 469.06 \u00b1 2.44 \u2502 19.78 \u00b1 0.45 \u2502
\u2502     Zjets     \u2502  -0.02 \u00b1 0.13  \u2502  2.6 \u00b1 2.22   \u2502 6.47 \u00b1 5.51  \u2502
\u2502    Zgamma     \u2502   0.0 \u00b1 0.0    \u2502   0.0 \u00b1 0.0   \u2502  0.3 \u00b1 0.29  \u2502
\u2502     ttbar     \u2502   0.0 \u00b1 0.0    \u2502  0.63 \u00b1 0.18  \u2502  0.28 \u00b1 0.1  \u2502
\u2502      WZ       \u2502   0.36 \u00b1 0.1   \u2502  1.79 \u00b1 0.23  \u2502 2.24 \u00b1 0.29  \u2502
\u2502      tZ       \u2502  0.01 \u00b1 0.01   \u2502  0.07 \u00b1 0.03  \u2502 0.06 \u00b1 0.02  \u2502
\u2502      ttZ      \u2502  1.24 \u00b1 0.08   \u2502  4.71 \u00b1 0.16  \u2502 5.74 \u00b1 0.18  \u2502
\u2502      tWZ      \u2502  0.57 \u00b1 0.11   \u2502  2.16 \u00b1 0.23  \u2502  2.5 \u00b1 0.24  \u2502
\u2502      VBS      \u2502  11.71 \u00b1 0.09  \u2502  6.4 \u00b1 0.08   \u2502 0.18 \u00b1 0.01  \u2502
\u2502      VH       \u2502  1.29 \u00b1 0.71   \u2502  5.76 \u00b1 1.4   \u2502 5.77 \u00b1 1.29  \u2502
\u2502    Others     \u2502  0.06 \u00b1 0.01   \u2502  0.4 \u00b1 0.13   \u2502 0.56 \u00b1 0.08  \u2502
\u251c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2524
\u2502   Bkg Tot.    \u2502 1790.82 \u00b1 5.94 \u2502 493.58 \u00b1 3.61 \u2502 43.89 \u00b1 5.7  \u2502
\u251c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u253c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2524
\u2502 Significance  \u2502   0.28 \u00b1 0.0   \u2502  0.42 \u00b1 0.0   \u2502  1.56 \u00b1 0.1  \u2502
\u2502 Combined Sig. \u2502   NaN \u00b1 0.0    \u2502  1.64 \u00b1 0.09  \u2502  NaN \u00b1 0.0   \u2502
\u2514\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2534\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2534\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2534\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2518
```
"""
print_sigtable(full_table) = pretty_table(
    full_table;
    header = ["", "SF-inZ", "SF-noZ", "DF"], 
    formatters = sigtable_fmt,
    body_hlines = [size(full_table, 1) - 3, size(full_table, 1) - 2],
    highlighters = hl_row([1], crayon"bold"), 
    alignment=:c,
)
