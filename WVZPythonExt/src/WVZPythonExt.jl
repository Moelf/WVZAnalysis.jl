module WVZPythonExt
using PythonCall, FHist, Serialization

export serial_to_root

function make_TH1D(h)
    np = pyimport("numpy")
    pyhist = pyimport("hist")

    bc = bincounts(h)
    hout = pyhist.Hist(pyhist.axis.Variable(np.array(binedges(h))), storage=pyhist.storage.Weight())
    va = h.sumw2
    hout[pybuiltins.Ellipsis] = np.stack([bc, va], axis=-1)
    return hout
end

"""
# Example
```
    julia> WVZPythonExt.find_CR_Njet_sys(hs)
    181-element Vector{String}:
     "JET_JER_EffectiveNP_5__1up"
     "MUON_EFF_ISO_STAT__1up"
     "JET_EtaIntercalibration_Modelling__1down"
     "FT_EFF_Eigen_B_0__1down"
     "JET_EtaIntercalibration_NonClosure_posEta__1up"
```
 """
function find_CR_Njet_sys(Hs)
    ks = string.(keys(Hs))
    filter!(ks) do s
        startswith(s, "ZZCR0j__Njet__")
    end
    return map(ks) do s
        rg = findfirst("Njet__", s)
        s[last(rg)+1:end]
    end
end

function combined_ZZCR(HS)
    ks = find_CR_Njet_sys(HS)
    Dict(
        Symbol("ZZCR__Njet__$k") =>
        (HS[Symbol("ZZCR0j__Njet__$k")] + HS[Symbol("ZZCR1j__Njet__$k")] + HS[Symbol("ZZCR2plusj__Njet__$k")])
        for k in ks
    )
end
function combined_ttZCR(HS)
    ks = find_CR_Njet_sys(HS)
    Dict(
        Symbol("ttZCR__Njet__$k") =>
        (HS[Symbol("ttZCR0j__Njet__$k")] + HS[Symbol("ttZCR1j__Njet__$k")] + HS[Symbol("ttZCR2plusj__Njet__$k")])
        for k in ks
    )
end

function serial_to_root(p)
    up = pyimport("uproot")
    isdir(p) || error("$p is not a directory")
    for fname in readdir(p)
        m = match(r"(.*)\.jlserialize", fname)
        isnothing(m) && continue
        tag = m[1]
        path = joinpath(p, "$(tag).jlserialize")
        !isfile(path) && continue
        Hs = deserialize(path)
        if tag == "ZZ"
            for ext in ("0j", "1j", "2plusj")
                pywith(up.recreate(joinpath(p, "$(tag)$(ext).root"))) do file
                    for (k,v) in Hs
                        k_str = string(k)
                        if !(contains(k_str, "CR") && contains(k_str, "Njet"))
                            # other histograms are already combined
                            # file[k_str] = make_TH1D(v)
                        elseif contains(k_str, ext)
                            file[replace(k_str, ext=>"")] = make_TH1D(v)
                        end
                    end
                end
            end
        end
        pywith(up.recreate(joinpath(p, "$(tag).root"))) do file
            for (k,v) in Hs
                k_str = string(k)
                if !(contains(k_str, "CR") && contains(k_str, "Njet"))
                    file[k_str] = make_TH1D(v)
                end
            end
            for (k,v) in combined_ZZCR(Hs)
                k_str = string(k)
                file[k_str] = make_TH1D(v)
            end
            for (k,v) in combined_ttZCR(Hs)
                k_str = string(k)
                file[k_str] = make_TH1D(v)
            end
        end
    end
end
end
