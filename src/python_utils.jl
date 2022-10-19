using PythonCall, Serialization

function make_TH1D(h)
    np = pyimport("numpy")
    pyhist = pyimport("hist")

    hout = pyhist.Hist(pyhist.axis.Regular(100, 0, 1), storage=pyhist.storage.Weight())
    bc = bincounts(h)
    va = h.sumw2
    hout[pybuiltins.Ellipsis] = np.stack([bc, va], axis=-1)
    return hout
end


function serial_to_root(p)
    up = pyimport("uproot")
    isdir(p) || error("$p is not a directory")
    for tag in WVZAnalysis.ALL_TAGS
        path = joinpath(p, "$(tag)_pmap.jlserialize")
        !isfile(path) && continue
        Hs = deserialize(path)
        pywith(up.recreate(joinpath(p, "$(tag).root"))) do file
            for (k,v) in Hs
                file[string(k)] = make_TH1D(v)
            end
        end
    end
end
