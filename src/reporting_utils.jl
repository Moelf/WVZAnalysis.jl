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

function significance_matrix(proc_names; recreate)
    M = mapreduce(vcat, proc_names) do tag
        ## re-make
        if recreate
            res = WVZAnalysis.sfsys(tag; NN_hist=true)
            serialize("/data/jiling/WVZ/v2.3_hists/$(tag).jlserialize", res)
        else
            #load from serialization
            res = deserialize("/data/jiling/WVZ/v2.3_hists/$(tag).jlserialize")
        end
        N = nbins(res[:DF__NN__NOMINAL])
        hists = rebin.([res[:SFinZ__NN__NOMINAL], res[:SFnoZ__NN__NOMINAL], res[:DF__NN__NOMINAL]], N)
        permutedims(hists) 
    end
    M
end

function significance_table(proc_names = ALL_TAGS; recreate=false)
    body = significance_matrix(proc_names; recreate)
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

    full_body = [collect(proc_names)
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
