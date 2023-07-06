module WVZXGBoostExt

using Dates, ROCCurves, Arrow, DataFrames, CSV
import XGBoost
using XGBoost: DMatrix

const ALL_TAGS = ["Signal", "ZZ", "Zjets", "ttbar", "WZ", "tZ", "ttZ", "tWZ", "VBS", "VH", "Others"];

const useful_features = [
    :leptonic_HT, :MET, :Zlep1_dphi,  :Wlep1_pt, :total_HT, :Zlep2_dphi, :Zlep2_eta, 
    :Njet, :Wlep2_eta,  :Zlep2_pt, :METSig, :other_mass, :Wlep1_dphi, 
    :Zlep1_pt,  :mass_4l, :pt_4l, :Zlep1_eta, :HT, 
    :Wlep1_eta, :Wlep2_dphi, :Zcand_mass, :Wlep2_pt, :MET_dPhi, :event, :SR,
];


function load_all_arrow(path) 
    df_all = mapreduce(vcat, ALL_TAGS) do tag
        df = DataFrame(Arrow.Table(joinpath(path, "$tag.arrow")))
        select!(df, useful_features)
        df.is_signal .= tag == "Signal" ? 1.f0 : 0.f0
        @. df.event = mod(df.event, 5)
        df
    end;

    return df_all
end

function train_and_log(df_all; output_dir, tree_method="gpu_hist")
    if !isdir(output_dir)
        @info "Creating dir $output_dir"
        mkdir(output_dir)
    end
    for (SR, SR_num) in zip([:SFinZ, :SFnoZ, :DF], 0:2)
        df_folds = [df_all[df_all.event .== i .&& df_all.SR .== SR_num, :] for i = 0:4]
            for fold = 0:4
                @info "Training $SR fold $fold"
                df_train = reduce(vcat, df_folds[Not(fold+1)])
                x_train = df_train[:, Not([:SR, :is_signal, :event])]
                y_train = df_train.is_signal;

                df_eval = df_folds[fold+1]
                x_eval = df_eval[:, Not([:SR, :is_signal, :event])]
                y_eval = df_eval.is_signal;

                res, log = xgboost_log(DMatrix(x_train, y_train);
                num_round=250,
                testdm = DMatrix(x_eval, y_eval),
                eta = 0.05,
                max_depth=8,
                colsample_bytree=0.8, 
                lambda = 1.0,
                alpha = 0.0,
                tree_method
                );

                # these two functions use different order convention, cursed if you ask me
                XGBoost.save(res, joinpath(output_dir, "$(SR)$(fold).model"))
                CSV.write(joinpath(output_dir, "$(SR)$(fold).log"), log)
            end
        end
    end



function xgboost_log(traindm::DMatrix;
                        testdm::Any=nothing,
                        num_round,
                        kw...
                    )
    
    Xy = XGBoost.DMatrix(traindm)
    b = XGBoost.Booster(Xy; kw...)
    update_feature_names::Bool=false
    if typeof(testdm)==DMatrix
        watchlist=Dict("train"=>traindm , "test"=>testdm)
    else
        watchlist=Dict("train"=>traindm)
    end
    names = collect(Iterators.map(string, keys(watchlist)))
    watch = collect(Iterators.map(x -> x.handle, values(watchlist)))
    thelog = Vector{String}(undef,0)
    for j in 1:num_round
        XGBoost.xgbcall(XGBoost.XGBoosterUpdateOneIter, b.handle, j, Xy.handle)
        o = Ref{Ptr{Int8}}()
        XGBoost.xgbcall(XGBoost.XGBoosterEvalOneIter, b.handle, j, watch, names, length(watch), o)
        push!(thelog,unsafe_string(o[]))
        XGBoost._maybe_update_feature_names!(b, Xy, update_feature_names)
    end
    return (booster=b , log=parsethelog(thelog))
end

function parsethelog(thelog::Vector{String})
    nr=length(thelog)
    neval= length(findall(":",thelog[1]))
    cstr=split(replace(replace(thelog[1],"\t"=>","),":"=>","),",")
    evalnames=Vector{String}(undef,0)
    for c in 1:neval
        push!(evalnames, cstr[2*c])
    end
    vals=zeros(nr,neval)
    rnd=zeros(Int,nr)
    for r in 1:nr
        l1=replace(thelog[r],"\t"=>",")
        l2=replace(l1,":"=>",")
        l3=split(l2,",")
        rnd[r]= parse(Int64,SubString(l3[1],2,length(l3[1])-1))
        for c in 1:neval
            vals[r,c]= parse(Float64,l3[1+2*c])
        end
    end
    valdf=hcat(DataFrame(iteration=rnd),DataFrame(vals, evalnames))
    return valdf
end

### for visualization, need CairoMakie
function loss_history(region)
    f = Figure()
    ax = Axis(f[1,1]; title=region, ylabel = "avg. log loss", xlabel = "N trees") #limits=(nothing, nothing, 0.1,0.5))
    for fold = 0:4
        df = CSV.read("/data/jiling/WVZ/v2.3-2023_06_12_hists/$(region)$(fold).log", DataFrame)
        fold+=1
        lines!(ax, df."train-rmse"; color = Cycled(fold), label = "Train fold: $fold")
        lines!(ax, df."test-rmse"; color = Cycled(fold), linestyle=:dash, label = "Test fold: $fold")
    end
    axislegend(; position=:rt, nbanks=2)
    f
end


end # module WVZXGBoostExt
