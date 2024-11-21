function populate_hist!(dict, shape_variation, symbols, bins, sfsys)
    for n in symbols
        dict[Symbol(n, :__, shape_variation)] = Hist1D(Float64; bins, overflow=true)
        !sfsys && continue
        for (_,vs) in SF_BRANCH_DICT
            for v in vs
                dict[Symbol(n, :__, v)] = Hist1D(Float64; bins, overflow=true)
            end
        end
    end
end

function init_BDT()
    SFinZs = Tuple(
              XGBoost.Booster(XGBoost.DMatrix[],  model_file = joinpath(BDT_MODEL_DIR, "SFinZ$i.model"), nthread=1)
              for i = 0:4
             )
    SFnoZs = Tuple(
              XGBoost.Booster(XGBoost.DMatrix[],  model_file = joinpath(BDT_MODEL_DIR, "SFnoZ$i.model"), nthread=1)
              for i = 0:4
             )
    DFs = Tuple(
              XGBoost.Booster(XGBoost.DMatrix[],  model_file = joinpath(BDT_MODEL_DIR, "DF$i.model"), nthread=1)
              for i = 0:4
             )
    return function f(ary; fold, region)
        bst = if region == :SFinZ
            SFinZs[fold]
        elseif region == :SFnoZ
            SFnoZs[fold]
        elseif region == :DF
            DFs[fold]
        end
        return Float32(predict(bst, permutedims(ary))[1])
    end
end


"""
    struct AnalysisTask
        path::String
        sumWeight::Float64
        isdata::Bool = false
        BDT_hist::Bool = true
        arrow_making::Bool = false
        sfsys::Bool = false
        shape_variation::String = "NOMINAL"
        controlregion::Symbol = :none
        require_VHSig::Union{Nothing, Bool} = nothing
        sel_Njet::Union{Nothing, Int} = nothing
    end

Fully define a task to be run on an executor, by calling `main_looper(task::AnalysisTask)`; see also
[`main_looper`](@ref).

You most likely don't need to construct it manually, see [`prep_tasks`](@ref):

# Example

```julia
julia> prep_tasks("Signal") |> first
path="/data/jiling/WVZ/v2.3/user.jiling.WVZ_v2.3sf.363507.e6379_s3126_r10201_p4434_ANALYSIS.root/user.jiling.29896106._000001.ANALYSIS.root"
sumWeight=13812.79638671875
isdata=false
BDT_hist=true
arrow_making=false
sfsys=false
shape_variation="NOMINAL"
controlregion=:none

julia> prep_tasks("Signal"; arrow_making=true) |> first
ERROR: can't do produce arrow and NN histograms at the same time
Stacktrace:
...
..

julia> prep_tasks("Signal"; arrow_making=true, BDT_hist=false) |> first
path="/data/jiling/WVZ/v2.3/user.jiling.WVZ_v2.3sf.363507.e6379_s3126_r10201_p4434_ANALYSIS.root/user.jiling.29896106._000001.ANALYSIS.root"
sumWeight=13812.79638671875
isdata=false
BDT_hist=false
arrow_making=true
sfsys=false
shape_variation="NOMINAL"
controlregion=:none
```

"""
struct AnalysisTask
    path::String
    sumWeight::Float64
    isdata::Bool
    BDT_hist::Bool
    arrow_making::Bool
    sfsys::Bool
    shape_variation::String
    controlregion::Symbol
    require_VHSig::Union{Nothing, Bool}
    sel_Njet::Union{Nothing, Int}
    function AnalysisTask(; path, sumWeight, isdata=false, BDT_hist=true,
            arrow_making=false, sfsys=false, shape_variation="NOMINAL",
            controlregion=:none, require_VHSig = nothing, sel_Njet = nothing)
        shapesys = (shape_variation != "NOMINAL")
        if sfsys && shapesys
            error("can't do SF systematic and Shape systematic at the same time")
        elseif arrow_making && BDT_hist
            error("can't produce arrow and NN histograms at the same time")
        else
            new(path, sumWeight, isdata, BDT_hist, arrow_making, sfsys, shape_variation, controlregion, require_VHSig, sel_Njet)
        end
    end
end

function Base.show(io::IO, a::AnalysisTask)
    for n in propertynames(a)
        print(io, "$n=")
        show(io, getproperty(a, n))
        println(io)
    end
end

function Base.filesize(a::AnalysisTask)
    filesize(a.path)
end



"""

# Example:

```julia
julia> dd = Dict(:Nlep => [], :lep1_pid=>[]);

julia> Nlep = 10;

julia> lep1_pid =11;

julia> @fill_dict! dd push! Nlep,lep1_pid;

# this expands to
push!(dd[:Nlep], Nlep)
push!(dd[:lep1_pid], lep1_pid)

julia> dd
Dict{Symbol, Vector{Any}} with 2 entries:
  :Nlep     => [10]
  :lep1_pid => [11]
```

"""
macro fill_dict!(dict, func, vars)
    vs = unique(vars.args)
    exs = [Expr(:call, func, Expr(:ref, dict, QuoteNode(v)), v) for v in vs]
    esc(Expr(:block, exs...))
end

macro fill_dict!(dict, wgt, func, vars)
    vs = unique(vars.args)
    exs = [Expr(:call, func, Expr(:ref, dict, QuoteNode(v)), v, wgt) for v in vs]
    esc(Expr(:block, exs...))
end

function extrac_dsid(str)
    return match(r"(\d{6})", str).captures[1] #dsid
end

"""
    root_dirs(tag::AbstractString; variation = "sf")

Reutrn a list of dir paths associated with a `tag` +sf or +shape
"""
function root_dirs(tag::AbstractString; variation = "sf")
    _LIST = joinpath(dirname(@__DIR__), "config/file_list.json") |> read |> JSON3.read
    folders = _LIST[tag]
    if lowercase(tag) == "data"
        sel1 = filter(readdir(MINITREE_DIR)) do folder_name
            contains(folder_name, "_13TeV") && contains(folder_name, "data")
        end
        return joinpath.(MINITREE_DIR, sel1)
    else
        dsids = unique(map(extrac_dsid, folders))
        # folder_name needs to contain any one element of the dsids
        sel2 = filter(readdir(MINITREE_DIR)) do folder_name
            any(occursin(folder_name), dsids) && contains(folder_name, variation)
        end
        return joinpath.(MINITREE_DIR, sel2)
    end
end

root_dirs(tags; variation = "sf") = mapreduce(x->root_dirs(x; variation), vcat, unique(tags))

"""
    Compute the sum of weights given path to a directory, and cache the result in a text file.
"""
function sumsumWeight(dir_path)
    cache = joinpath(dir_path, "sumsumWeight.txt")
    if isfile(cache)

        res = parse(Float64, readchomp(cache))
        return res
    end
    
    res = 0.0
    for f in readdir(dir_path; join=true)
        r = ROOTFile(f)
        if !haskey(r, "sumWeight")
            # in case of data
            res = 1.0
            break
        end
        res += r["sumWeight"][:fN][3]
    end
    # 2.4.1 has double counting, need divide by 2
    if contains(dir_path, "2.4.1sf")
        res /= 2
    end
    try
    open(cache, "w") do io
        println(io, res)
    end
    catch
    end
    return res
end


"""
    prep_tasks(tag; shape_variation="NOMINAL", scouting=false, kw...)
    
For construction a collection of `AnalysisTask`s given tag and job options,
check out [@AnalysisTask].

All posible tags can be found in `config/file_list.json`, there's also a convinient
variable `WVZAnalysis.ALL_TAGS` that keeps track of all processes in an exclusive mannar:

# Example

```julia
julia> WVZAnalysis.ALL_TAGS
("Signal", "ZZ", "Zjets", "Zgamma", "ttbar", "WZ", "tZ", "ttZ", "tWZ", "VBS", "VH", "Others")
```

```julia
julia> all_nominal_tasks = mapreduce(prep_tasks, vcat, WVZAnalysis.ALL_TAGS);

julia> length(all_nominal_tasks)
768
```
"""
function prep_tasks(tag; shape_variation="NOMINAL", scouting=false, kw...)
    sel_Njet = nothing
    if startswith(tag, "ZZ")
        sel_Njet = if startswith(tag, "ZZ0")
            0
        elseif startswith(tag, "ZZ1")
            1
        elseif startswith(tag, "ZZ2")
            2
        end
        tag="ZZ"
    end

    dirs = if shape_variation == "NOMINAL"
        root_dirs(tag; variation = "sf")
    else
        root_dirs(tag; variation = "shape")
    end

    isdata = (lowercase(tag) == "data")

    if scouting
        @info "scouting"
        dirs = first(dirs, 1)
    end

    require_VHSig = if contains(tag, "H_Signal")
        true
    elseif contains(tag, "H_Bkg")
        false
    else
        nothing
    end
    


    files = mapreduce(vcat, dirs) do d
        PATHS = dir_to_paths(d; scouting)
        sumWeight = sumsumWeight(d)
        if occursin(r"341450", d)
            sumWeight *= (11.007105/2.412326)
        end
        if occursin(r"341452", d)
            sumWeight *= (11.186538/2.409899)
        end
        if occursin(r"341454", d)
            sumWeight *= (11.101932/2.411167)
        end
        if occursin(r"345066", d)
            sumWeight *= (57.429000/3.368E-02)
        end
        [AnalysisTask(; path, sumWeight, isdata, shape_variation, require_VHSig, sel_Njet, kw...) for path in PATHS]
    end
    files
end

function dir_to_paths(dir_path; scouting = false)
    paths = filter!(endswith(".root"), readdir(dir_path; join = true))
    if scouting
        @info "scouting"
        paths = first(paths, 2)
    end
    return paths
end

"""
    arrow_making(tasks; mapper = map)

Take a collection of tasks, run them via `map` and `mergewith(append!)`.
Returns a `dict` of vectors representing the datas after filtering.
"""
function arrow_making(tasks; mapper = map)
    p = Progress(length(tasks))
    res = mapper(tasks) do t
        x = main_looper(t)
        next!(p)
        x
    end
    
    return reduce(mergewith(append!), res)
end


"""
    hist_main(tag; mapfun=robust_pmap, no_shape = false, output_dir, kw...)

Running the main looper for a tag (e.g. VH, ttZ) to produce histogram
kw... takes anything that AnalysisTask takes.
"""
function hist_main(tag; mapfun=robust_pmap, no_shape = false, output_dir, kw...)
    p = output_dir
    if !isdir(p)
        mkdir(p)
    end

    all_tasks = prep_tasks(tag; sfsys=true)
    if tag != "Data"
        if no_shape || tag in ("WH_Bkg", "ZH_Bkg")
            @info "-------------- $tag SF ------------ "
        else
            @info "-------------- $tag SF + shapes ------------ "
            shape_tasks = mapreduce(shape_variation -> prep_tasks(tag; shape_variation, sfsys=false), vcat,
                                    SHAPE_TREE_NAMES)
            sort!(shape_tasks; by = x->x.path)
            append!(all_tasks, shape_tasks)
        end
    else
        @info "-------------- !!! $tag !!! ------------ "
    end
    n = length(all_tasks)
    println("$n tasks in total")
    pm = Progress(n; dt=1.0)
    all_list = mapfun(all_tasks) do task
        _m = main_looper(task)
        next!(pm)
        return _m
    end
    # TODO
    # process_theory_syst!.(all_list)
    Hs = reduce(mergewith(+), all_list)

    serialize(joinpath(p, "$(tag).jlserialize"), Hs)
    Hs
end

"""
    arrow_main(tag; mapfun=OhMyThreads.tmap, output_dir, kw...)

Running the main looper for a tag (e.g. VH, ttZ) to produce arrow
kw... takes anything that AnalysisTask takes.

"""
function arrow_main(tag; mapfun=OhMyThreads.tmap, output_dir, kw...)
    p = output_dir
    if !isdir(p)
        mkdir(p)
    end

    all_tasks = prep_tasks(tag; arrow_making=true, BDT_hist=false)
    if tag != "Data"
        @info "-------------- $tag SF ------------ "
    else
        @info "-------------- !!! $tag !!! ------------ "
    end
    n = length(all_tasks)
    println("$n tasks in total")
    pm = Progress(n; dt=1.0)
    all_list = mapfun(all_tasks) do task
        _m = main_looper(task)
        next!(pm)
        return _m
    end

    Hs = reduce(mergewith(vcat), all_list)

    Arrow.write(joinpath(p, "$(tag).arrow"), Hs)
end

"""
    make_sfsys_wgt!(evt, wgt, nominal_branch, wgt_name, wgt_idx=1)

1. use `wgt_name` to find alternative branch names
2. For every variation, make a new entry in `wgt::Dict`.
3. Also multiply nominal for `wgt[:NOMINAL]`

"""
function make_sfsys_wgt!(evt, wgt, wgt_name, wgt_idx=1; pre_mask=Colon(), sfsys)
    orig_wgt = wgt[:NOMINAL]
    nominal_branch = getproperty(evt, wgt_name)
    nominal_new = reduce(*, @views(nominal_branch[pre_mask][wgt_idx]))
    if sfsys
        variation_names = SF_BRANCH_DICT[wgt_name]
        for key_name in variation_names
            full_name = Symbol(wgt_name, :__, key_name)
            var_branch = getproperty(evt, full_name)
            var_wgt = reduce(*, @views(var_branch[pre_mask][wgt_idx]))
            # default get orig_wgt when first access
            wgt[key_name] = get(wgt, key_name, orig_wgt) * var_wgt
        end
        for key_name in keys(wgt)
            key_name ∈ variation_names && continue
            wgt[key_name] *= nominal_new
        end
    else
        wgt[:NOMINAL] *= nominal_new
    end
    nothing
end

function cutflow_total!(cutflow_ptr, dict, wgt_dict; BDT_hist, shape_variation)
    if BDT_hist && shape_variation == "NOMINAL"
        cutflow_ptr[] += 1
        push!(dict[:CutFlow], cutflow_ptr[])
        push!(dict[:CutFlowWgt], cutflow_ptr[], wgt_dict[:NOMINAL])
    end
end
function cutflow_SRs!(cutflow_ptr, dict, wgt_dict; BDT_hist, SR, shape_variation)
    if BDT_hist && shape_variation == "NOMINAL"
        cutflow_ptr[] += 1
        if SR == 0
            push!(dict[:SFinZCutFlow], cutflow_ptr[])
            push!(dict[:SFinZCutFlowWgt], cutflow_ptr[], wgt_dict[:NOMINAL])
        elseif SR == 1
            push!(dict[:SFnoZCutFlow], cutflow_ptr[])
            push!(dict[:SFnoZCutFlowWgt], cutflow_ptr[], wgt_dict[:NOMINAL])
        else
            push!(dict[:DFCutFlow], cutflow_ptr[])
            push!(dict[:DFCutFlowWgt], cutflow_ptr[], wgt_dict[:NOMINAL])
        end
    end
end

function BDT_hist_init(; sfsys, shape_variation)
    if sfsys && (shape_variation != "NOMINAL")
        error("can't do sf systematics and shape systematics at the same time")
    end
    _dict = Dict{Symbol, Hist1D}()

    bins = 0:0.01:1
    populate_hist!(_dict, shape_variation, (:SFinZ__BDT, :SFnoZ__BDT, :DF__BDT), bins, sfsys)

    bins = 0:5:300
    populate_hist!(_dict, shape_variation, (:SFinZ__MET, :SFnoZ__MET, :DF__MET), bins, sfsys)

    bins = 0:5
    populate_hist!(_dict, shape_variation, (:SFinZ__Njet, :SFnoZ__Njet, :DF__Njet), bins, sfsys)
    populate_hist!(_dict, shape_variation, (:ZZCR__Njet, :ttZCR__Njet), bins, sfsys)

    bins = 0:10:1000
    populate_hist!(_dict, shape_variation, (:ttZCR__m4l, ), bins, sfsys)

    leptonic_HT_names = (Symbol(x, :__leptonic_HT) for x in (:SFinZ, :SFnoZ, :DF))
    HT_names = (Symbol(x, :__HT) for x in (:SFinZ, :SFnoZ, :DF))
    total_HT_names = (Symbol(x, :__total_HT) for x in (:SFinZ, :SFnoZ, :DF))
    MET_names = (Symbol(x, :__MET) for x in (:SFinZ, :SFnoZ, :DF))
    METsig_names = (Symbol(x, :__METsig) for x in (:SFinZ, :SFnoZ, :DF))
    METdphi_names = (Symbol(x, :__METdphi) for x in (:SFinZ, :SFnoZ, :DF))
    Z_mass_names = (Symbol(x, :__Z_mass) for x in (:SFinZ, :SFnoZ, :DF))
    WW_mass_names = (Symbol(x, :__WW_mass) for x in (:SFinZ, :SFnoZ, :DF))
    pt4l_names = (Symbol(x, :__pt4l) for x in (:SFinZ, :SFnoZ, :DF))
    Wl_pt1_names = (Symbol(x, :__Wl_pt1) for x in (:SFinZ, :SFnoZ, :DF))
    Wl_pt2_names = (Symbol(x, :__Wl_pt2) for x in (:SFinZ, :SFnoZ, :DF))
    Wl_eta1_names = (Symbol(x, :__Wl_eta1) for x in (:SFinZ, :SFnoZ, :DF))
    Wl_eta2_names = (Symbol(x, :__Wl_eta2) for x in (:SFinZ, :SFnoZ, :DF))
    Zl_pt1_names = (Symbol(x, :__Zl_pt1) for x in (:SFinZ, :SFnoZ, :DF))
    Zl_pt2_names = (Symbol(x, :__Zl_pt2) for x in (:SFinZ, :SFnoZ, :DF))
    Zl_eta1_names = (Symbol(x, :__Zl_eta1) for x in (:SFinZ, :SFnoZ, :DF))
    Zl_eta2_names = (Symbol(x, :__Zl_eta2) for x in (:SFinZ, :SFnoZ, :DF))
    Zl_dphi1_names = (Symbol(x, :__Zl_dphi1) for x in (:SFinZ, :SFnoZ, :DF))
    Zl_dphi2_names = (Symbol(x, :__Zl_dphi2) for x in (:SFinZ, :SFnoZ, :DF))
    Wl_dphi1_names = (Symbol(x, :__Wl_dphi1) for x in (:SFinZ, :SFnoZ, :DF))
    Wl_dphi2_names = (Symbol(x, :__Wl_dphi2) for x in (:SFinZ, :SFnoZ, :DF))


    bins = range(length=101, start=0, stop=1000);
    populate_hist!(_dict, shape_variation, leptonic_HT_names, bins, sfsys)
    populate_hist!(_dict, shape_variation, HT_names, bins, sfsys)
    populate_hist!(_dict, shape_variation, total_HT_names, bins, sfsys)
    bins = range(length=101, start=0, stop=400);
    populate_hist!(_dict, shape_variation, MET_names, bins, sfsys)
    bins = range(length=101, start=0, stop=100);
    populate_hist!(_dict, shape_variation, METsig_names, bins, sfsys)
    bins = range(length=101, start=-3.14, stop=3.14);
    populate_hist!(_dict, shape_variation, METdphi_names, bins, sfsys)
    bins = range(length=101, start=0, stop=400);
    populate_hist!(_dict, shape_variation, Z_mass_names, bins, sfsys)
    populate_hist!(_dict, shape_variation, WW_mass_names, bins, sfsys)

    bins = range(length=101, start=0, stop=1000);
    populate_hist!(_dict, shape_variation, pt4l_names, bins, sfsys)

    bins = range(length=101, start=0, stop=400);
    populate_hist!(_dict, shape_variation, Wl_pt1_names, bins, sfsys)
    populate_hist!(_dict, shape_variation, Wl_pt2_names, bins, sfsys)

    bins = range(length=31, start=-3, stop=3);
    populate_hist!(_dict, shape_variation, Wl_eta1_names, bins, sfsys)
    populate_hist!(_dict, shape_variation, Wl_eta2_names, bins, sfsys)

    bins = range(length=101, start=0, stop=400);
    populate_hist!(_dict, shape_variation, Zl_pt1_names, bins, sfsys)
    populate_hist!(_dict, shape_variation, Zl_pt2_names, bins, sfsys)

    bins = range(length=31, start=-3, stop=3);
    populate_hist!(_dict, shape_variation, Zl_eta1_names, bins, sfsys)
    populate_hist!(_dict, shape_variation, Zl_eta2_names, bins, sfsys)

    bins = range(length=41, start=-3.14, stop=3.14);
    populate_hist!(_dict, shape_variation, Zl_dphi1_names, bins, sfsys)
    populate_hist!(_dict, shape_variation, Zl_dphi2_names, bins, sfsys)
    populate_hist!(_dict, shape_variation, Wl_dphi1_names, bins, sfsys)
    populate_hist!(_dict, shape_variation, Wl_dphi2_names, bins, sfsys)

    _dict[:CutFlow] = Hist1D(Int; bins=1:20)
    _dict[:CutFlowWgt] = Hist1D(Float64; bins=1:20)
    _dict[:SFinZCutFlow] = Hist1D(Int; bins=1:20)
    _dict[:SFinZCutFlowWgt] = Hist1D(Float64; bins=1:20)
    _dict[:SFnoZCutFlow] = Hist1D(Int; bins=1:20)
    _dict[:SFnoZCutFlowWgt] = Hist1D(Float64; bins=1:20)
    _dict[:DFCutFlow] = Hist1D(Int; bins=1:20)
    _dict[:DFCutFlowWgt] = Hist1D(Float64; bins=1:20)

    return _dict
end

function arrow_init()
    return Dict(
                           #lep 1,2,3,4 ordered by decending Pt 
                           :SR => Int32[],
                           :event => UInt64[],
                           :sr_SF_inZ => Int32[],
                           :sr_SF_noZ => Int32[],
                           :sr_DF => Int32[],
                           :cr_ZZ => Int32[],
                           :cr_ttZ => Int32[],
                           :Nlep => Int32[],
                           :lep1_pid => Int32[], 
                           :lep2_pid => Int32[],
                           :lep3_pid => Int32[],
                           :lep4_pid => Int32[],
                           :pt_1 => Float32[],
                           :pt_2 => Float32[],
                           :pt_3 => Float32[],
                           :pt_4 => Float32[],
                           :eta_1 => Float32[],
                           :eta_2 => Float32[],
                           :eta_3 => Float32[],
                           :eta_4 => Float32[],
                           :phi_1 => Float32[],
                           :phi_2 => Float32[],
                           :phi_3 => Float32[],
                           :phi_4 => Float32[],
                           :Njet => Int32[],
                           :mass_4l => Float32[],
                           :Zcand_mass => Float32[],
                           :other_mass => Float32[],
                           :MET => Float32[],
                           :METSig => Float32[],
                           :METPhi => Float32[],
                           :MET_dPhi => Float32[],
                           :leptonic_HT => Float32[],
                           :HT => Float32[],
                           :total_HT => Float32[],
                           :Zlep1_pt => Float32[],
                           :Zlep1_eta => Float32[],
                           :Zlep1_phi => Float32[],
                           :Zlep1_dphi => Float32[],
                           :Zlep1_pid => Int32[],
                           :Zlep2_pt => Float32[],
                           :Zlep2_eta => Float32[],
                           :Zlep2_phi => Float32[],
                           :Zlep2_dphi => Float32[],
                           :Zlep2_pid => Int32[],
                           :Wlep1_pt => Float32[],
                           :Wlep1_eta => Float32[],
                           :Wlep1_phi => Float32[],
                           :Wlep1_dphi => Float32[],
                           :Wlep1_pid => Int32[],
                           :Wlep2_pt => Float32[],
                           :Wlep2_eta => Float32[],
                           :Wlep2_phi => Float32[],
                           :Wlep2_dphi => Float32[],
                           :Wleps_deta => Float32[],
                           :Wlep2_pid => Int32[],
                           :pt_4l => Float32[],
                           :wgt => Float64[],
                           :mcGenWgt => Float32[],
                           :jet_pt_1 => Float32[],
                           :jet_pt_2 => Float32[],
                           :jet_pt_3 => Float32[],
                           :jet_pt_4 => Float32[],
                           :jet_eta_1 => Float32[],
                           :jet_eta_2 => Float32[],
                           :jet_eta_3 => Float32[],
                           :jet_eta_4 => Float32[],
                           :jet_phi_1 => Float32[],
                           :jet_phi_2 => Float32[],
                           :jet_phi_3 => Float32[],
                           :jet_phi_4 => Float32[],
                           :jet_m_1 => Float32[],
                           :jet_m_2 => Float32[],
                           :jet_m_3 => Float32[],
                           :jet_m_4 => Float32[],
                           :jet_m_4 => Float32[],
                           :v_j_btagCont => Vector{Int32}[],
                           :v_j_btag60 => Vector{Bool}[],
                           :v_j_btag70 => Vector{Bool}[],
                           :v_j_btag77 => Vector{Bool}[],
                           :v_j_btag85 => Vector{Bool}[],
                           :jet_btagCont_1 => Int32[],
                           :jet_btagCont_2 => Int32[],
                           :jet_btagCont_3 => Int32[],
                           :jet_btagCont_4 => Int32[],
                          )
end
