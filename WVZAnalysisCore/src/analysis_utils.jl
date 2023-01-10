
"""
extract dsid from a file name, used to match with systematic files
"""
const MINITREE_DIR = Ref("/data/jiling/WVZ/v2.3")
const ANALYSIS_DIR = Ref("/data/jiling/WVZ/v2.3_hists")

const ONNX_MODEL_PATH = Ref("/data/grabanal/NN/NN_08_23.onnx")
const BDT_MODEL_PATH = Ref(joinpath(dirname(@__DIR__), "BDT_models"))

function init_ONNX()
    model=ONNX.load(ONNX_MODEL_PATH[], zeros(Float32, 30, 1))
    rescaling_parameters = joinpath(dirname(@__DIR__), "config/NN_08_23_rescaling_parameters.json") |> read |> JSON3.read
    NN_order = ("HT", "MET", "METPhi", "METSig", "Njet", "Wlep1_dphi", "Wlep1_eta",
                "Wlep1_phi", "Wlep1_pt", "Wlep2_dphi", "Wlep2_eta", "Wlep2_phi",
                "Wlep2_pt", "Zcand_mass", "Zlep1_dphi", "Zlep1_eta", "Zlep1_phi",
                "Zlep1_pt", "Zlep2_dphi", "Zlep2_eta", "Zlep2_phi", "Zlep2_pt",
                "leptonic_HT", "mass_4l", "other_mass", "pt_4l", "total_HT",
                "sr_SF_inZ", "sr_SF_noZ", "sr_DF")
    return model, 
    [rescaling_parameters["scale"][name] for name in NN_order],
    [rescaling_parameters["min"][name] for name in NN_order]
end

function init_BDT()
    SFinZs = Tuple(
              XGBoost.Booster(XGBoost.DMatrix[],  model_file = joinpath(BDT_MODEL_PATH[], "SFinZ$i.model"))
              for i = 0:4
             )
    SFnoZs = Tuple(
              XGBoost.Booster(XGBoost.DMatrix[],  model_file = joinpath(BDT_MODEL_PATH[], "SFnoZ$i.model"))
              for i = 0:4
             )
    DFs = Tuple(
              XGBoost.Booster(XGBoost.DMatrix[],  model_file = joinpath(BDT_MODEL_PATH[], "DF$i.model"))
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
        NN_hist::Bool = true
        arrow_making::Bool = false
        sfsys::Bool = false
        shape_variation::String = "NOMINAL"
        controlregion::Symbol = :none
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
NN_hist=true
arrow_making=false
sfsys=false
shape_variation="NOMINAL"
controlregion=:none

julia> prep_tasks("Signal"; arrow_making=true) |> first
ERROR: can't do produce arrow and NN histograms at the same time
Stacktrace:
...
..

julia> prep_tasks("Signal"; arrow_making=true, NN_hist=false) |> first
path="/data/jiling/WVZ/v2.3/user.jiling.WVZ_v2.3sf.363507.e6379_s3126_r10201_p4434_ANALYSIS.root/user.jiling.29896106._000001.ANALYSIS.root"
sumWeight=13812.79638671875
isdata=false
NN_hist=false
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
    NN_hist::Bool
    arrow_making::Bool
    sfsys::Bool
    shape_variation::String
    controlregion::Symbol
    function AnalysisTask(; path, sumWeight, isdata=false, NN_hist=true,
            arrow_making=false, sfsys=false, shape_variation="NOMINAL",
            controlregion=:none)
        shapesys = (shape_variation != "NOMINAL")
        if sfsys && shapesys
            error("can't do SF systematic and Shape systematic at the same time")
        elseif arrow_making && NN_hist
            error("can't produce arrow and NN histograms at the same time")
        else
            new(path, sumWeight, isdata, NN_hist, arrow_making, sfsys, shape_variation, controlregion)
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


function NN_calc(model, scales, minimums, NN_input)::Float32
    @. NN_input = fma(NN_input, scales, minimums)
    return Ghost.play!(model, NN_input)[1]
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
        sel1 = filter(readdir(MINITREE_DIR[])) do folder_name
            contains(folder_name, "_13TeV") && contains(folder_name, "data")
        end
        return joinpath.(MINITREE_DIR[], sel1)
    else
        dsids = unique(map(extrac_dsid, folders))
        # folder_name needs to contain any one element of the dsids
        sel2 = filter(readdir(MINITREE_DIR[])) do folder_name
            any(occursin(folder_name), dsids) && contains(folder_name, variation)
        end
        return joinpath.(MINITREE_DIR[], sel2)
    end
end

"""
multiple tags
"""
root_dirs(tags; variation = "sf") = mapreduce(x->root_dirs(x; variation), vcat, unique(tags))

function sumsumWeight(dir_path)::Float64
    cache = joinpath(dir_path, "sumsumWeight.txt")
    if isfile(cache)
        return parse(Float64, readchomp(cache))
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
    open(cache, "w") do io
        println(io, res)
    end
    return res
end

function _runwork(tasks; mapper=map)
    println("processing $(length(tasks)) root files in total.")
    P = Progress(length(tasks))
    s = mapper(tasks) do t
        main_looper(t)
        next!(P)
    end
    return reduce(mergewith(+), s)
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
    dirs = if shape_variation == "NOMINAL"
        root_dirs(tag; variation = "sf")
    else
        root_dirs(tag; variation = "shape")
    end

    isdata = (lowercase(tag) == "data")

    if scouting
        @info "scounting"
        dirs = first(dirs, 1)
    end

    files = mapreduce(vcat, dirs) do d
        PATHS = dir_to_paths(d; scouting)
        sumWeight = sumsumWeight(d)
        if occursin(r"346645|346646|346647", d)
            @show d
            sumWeight *= 2.745e-4
        end
        [AnalysisTask(; path, sumWeight, isdata, shape_variation, kw...) for path in PATHS]
    end
    files
end

function sfsys(tag; scouting=false, kw...)
    isdata = (lowercase(tag) == "data")
    tasks = prep_tasks(tag)
    return _runwork(tasks; isdata, kw...)
end

function sfsys_pmap(tag)
    isdata = (lowercase(tag) == "data")
    tasks = prep_tasks(tag)
    s = @showprogress @distributed mergewith((.+)) for t in tasks
        main_looper(t)
    end
    s
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
    arrow_making(tasks)

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
    hist_root(tag; mapper=pmap, no_shape=false, output_dir, kw...)

Similar to the one without `pmap`, except uses pmap for everything. Checkout `ClusterManager.jl`
and be sure to have a handful of workers before running this.
"""
function hist_root(tag; mapper=pmap, no_shape = false, output_dir, kw...)
    p = output_dir
    if !isdir(p)
        mkdir(p)
    end
    @show nworkers()

    all_tasks = prep_tasks(tag; sfsys=true)
    if tag != "Data"
        if no_shape
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
    println("$(length(all_tasks)) tasks in total")
    prog = Progress(length(all_tasks))
    all_list = mapper(all_tasks) do task
        _m = main_looper(task)
        next!(prog)
        return _m
    end
    # TODO
    # process_theory_syst!.(all_list)
    Hs = reduce(mergewith(+), all_list)

    serialize(joinpath(p, "$(tag).jlserialize"), Hs)
    Hs
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
