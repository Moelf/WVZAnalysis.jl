function NN_hist_init(; sfsys, shape_variation)
    if sfsys && (shape_variation != "NOMINAL")
        error("can't do sf systematics and shape systematics at the same time")
    end
    _dict = Dict{Symbol, Hist1D}()
    bins = 0:0.1:1
    for n in (:SFinZ__NN, :SFnoZ__NN, :DF__NN)
        _dict[Symbol(n, :__, shape_variation)] = Hist1D(Float64; bins, overflow=true)

        !sfsys && continue
        for (_,vs) in SF_BRANCH_DICT
            for v in vs
                _dict[Symbol(n, :__, v)] = Hist1D(Float64; bins, overflow=true)
            end
        end
    end
    bins = 0:5:300
    for n in (:SFinZ__MET, :SFnoZ__MET, :DF__MET)
        _dict[Symbol(n, :__, shape_variation)] = Hist1D(Float64; bins, overflow=true)

        !sfsys && continue
        for (_,vs) in SF_BRANCH_DICT
            for v in vs
                _dict[Symbol(n, :__, v)] = Hist1D(Float64; bins, overflow=true)
            end
        end
    end
    for n in (:CR__yield,)
        _dict[Symbol(n, :__, shape_variation)] = Hist1D(Float64; bins=0:1, overflow=true)

        !sfsys && continue
        for (_,vs) in SF_BRANCH_DICT
            for v in vs
                _dict[Symbol(n, :__, v)] = Hist1D(Float64; bins=0:1, overflow=true)
            end
        end
    end

    return _dict
end

function kinematic_hist_init()
    _dict = Dict(
                 :pt_1 => Hist1D(Float64; bins=0:5:300, overflow=true),
                 :pt_2 => Hist1D(Float64; bins=0:5:300, overflow=true),
                 :pt_3 => Hist1D(Float64; bins=0:5:300, overflow=true),
                 :pt_4 => Hist1D(Float64; bins=0:5:300, overflow=true),
                 :eta_1 => Hist1D(Float64; bins=-3:0.2:3, overflow=true),
                 :eta_2 => Hist1D(Float64; bins=-3:0.2:3, overflow=true),
                 :eta_3 => Hist1D(Float64; bins=-3:0.2:3, overflow=true),
                 :eta_4 => Hist1D(Float64; bins=-3:0.2:3, overflow=true),
                 :mass_4l => Hist1D(Float64; bins=0:30:1000, overflow=true),
                 :Zcand_mass => Hist1D(Float64; bins=50:2:150, overflow=true),
                 :other_mass => Hist1D(Float64; bins=0:4:500, overflow=true),
                 :MET => Hist1D(Float64; bins=bins=0:1:1000, overflow=true),
                 :METSig => Hist1D(Float64; bins=0:0.2:20, overflow=true),
                 :HT => Hist1D(Float64; bins=0:10:800, overflow=true),
                 :leptonic_HT => Hist1D(Float64; bins=0:10:800, overflow=true),
                 :total_HT => Hist1D(Float64; bins=0:10:2000, overflow=true),
                 :SR => Hist1D(Float64; bins=0:3),
                 :Z_eta => Hist1D(Float64; bins=-4:0.4:4, overflow=true),
                 :Z_phi => Hist1D(Float64; bins=-4:0.4:4, overflow=true),
                 :Z_pt => Hist1D(Float64; bins=0:8:600, overflow=true),
                 :Z_rapidity => Hist1D(Float64; bins=-4:0.4:4, overflow=true),
                 :total_events => Hist1D(Float64; bins=0:1, overflow=true),
                 :NN_score => Hist1D(Float64; bins=0:0.01:1, overflow=true),
                 :Njet => Hist1D(Float64; bins=0:10,overflow=true),
                )

    return _dict
end
