using PrettyTables, Serialization, Measurements

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
        if recreate
            res = WVZAnalysis.sfsys(tag; NN_hist=true)
            serialize("/data/jiling/WVZ/v2.3_hists/$(tag).jlserialize", res)
        else
            #load from serialization
            res = deserialize("/data/jiling/WVZ/v2.3_hists/$(tag).jlserialize")
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

    ## Example

```julia-repl
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
print_sigtable(full_table) = pretty_table(
    full_table;
    header = ["", "SF-inZ", "SF-noZ", "DF"], 
    formatters = sigtable_fmt,
    body_hlines = [size(full_table, 1) - 3, size(full_table, 1) - 2],
    highlighters = hl_row([1], crayon"bold"), 
    alignment=:c,
)
