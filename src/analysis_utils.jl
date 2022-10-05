
"""
extract dsid from a file name, used to match with systematic files
"""
const MINITREE_DIR = Ref("/data/jiling/WVZ/v2.3")

# const ONNX_MODEL_PATH = Ref("/data/grabanal/NN/NN_08_23.onnx")
const BDT_MODEL_PATH = Ref("/data/jiling/WVZ/v2.3-beta2_arrow/xgb_2022-09-27.model")

# function init_ONNX()
#     model=ONNX.load(ONNX_MODEL_PATH[], zeros(Float32, 30, 1))
#     rescaling_parameters = joinpath(dirname(@__DIR__), "config/NN_08_23_rescaling_parameters.json") |> read |> JSON3.read
#     NN_order = ("HT", "MET", "METPhi", "METSig", "Njet", "Wlep1_dphi", "Wlep1_eta",
#                 "Wlep1_phi", "Wlep1_pt", "Wlep2_dphi", "Wlep2_eta", "Wlep2_phi",
#                 "Wlep2_pt", "Zcand_mass", "Zlep1_dphi", "Zlep1_eta", "Zlep1_phi",
#                 "Zlep1_pt", "Zlep2_dphi", "Zlep2_eta", "Zlep2_phi", "Zlep2_pt",
#                 "leptonic_HT", "mass_4l", "other_mass", "pt_4l", "total_HT",
#                 "sr_SF_inZ", "sr_SF_noZ", "sr_DF")
#     return model, 
#     [rescaling_parameters["scale"][name] for name in NN_order],
#     [rescaling_parameters["min"][name] for name in NN_order]
# end
#

"""
    Base.@kwdef struct AnalysisTask
        path::String
        sumWeight::Float64
        isdata::Bool = false
        NN_hist::Bool = true
        arrow_making::Bool = false
        sfsys::Bool = false
        shape_variation::String = "NOMINAL"
        controlregion::Symbol = :none
    end

Fully define a task to be run on an executor, via `main_looper(job::AnalysisTask)`
"""
Base.@kwdef struct AnalysisTask
    path::String
    sumWeight::Float64
    isdata::Bool = false
    NN_hist::Bool = true
    arrow_making::Bool = false
    sfsys::Bool = false
    shape_variation::String = "NOMINAL"
    controlregion::Symbol = :none
end

function init_BDT()
    bst = Booster(; model_file = BDT_MODEL_PATH[])
    return ary->Float32(predict(bst, permutedims(ary))[1])
end


# function NN_calc(model, scales, minimums, NN_input)::Float32
#     @inbounds @simd for i in eachindex(NN_input, scales, minimums)
#         NN_input[i] = fma(NN_input[i], scales[i], minimums[i])
#     end
#     return Ghost.play!(model, NN_input)[1]
# end


"""

Example:

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
    s = @showprogress map(tasks) do t
        main_looper(t)
    end
    return mergewith((.+), s...)
end

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
        @info "scounting"
        paths = first(files, 2)
    end
    return paths
end

function arrow_making(tag)
    dirs = root_dirs(tag; variation = "sf")
    isdata = lowercase(tag) == "data"
    prog = Progress(mapreduce(length∘readdir, +, dirs), 0.5)
    println("$tag starting:")
    res = ThreadsX.map(dirs) do d
        arrow_making_dir(d; prog, isdata)
    end
    foldl((x,y) -> mergewith(append!, x, y), res)
    first(res)
end

function arrow_making_dir(dir_path; prog = nothing, isdata, kw...)
    files = filter!(endswith(".root"), readdir(dir_path; join = true))
    sumWeight = sumsumWeight(files)
    res = map(files) do F
        r = WVZAnalysis.main_looper(F, sumWeight; arrow_making=true, isdata, kw...)
        if prog !== nothing 
            next!(prog)
        end
        r
    end
    foldl((x,y) -> append!.(x,y), res)
    first(res)
end


function hist_root(tag; kw...)
    p = "/data/jiling/WVZ/v2.3_hists_uproot_oct4"
    if !isdir(p)
        mkdir(p)
    end
    sf1 = sfsys(tag; NN_hist=true, kw...);
    shape_list = [WVZAnalysis.shapesys(tag, variation; NN_hist=true, kw...)
                  for variation in WVZAnalysis.SHAPE_TREE_NAMES
            ]
    Hs = merge(sf1, shape_list...)
    serialize(joinpath(p, "$(tag).jlserialize"), Hs)
    Hs
end

function hist_root_pmap(tag; kw...)
    p = "/data/jiling/WVZ/v2.3_hists_uproot_oct5"
    if !isdir(p)
        mkdir(p)
    end
    @show nworkers()
    @info "-------------- $tag SF begin ------------ "
    sf_tasks = WVZAnalysis.prep_tasks(tag)
    println("$(length(sf_tasks)) tasks in total")
    sf_list = @showprogress pmap(main_looper, sf_tasks)
    sf_hist = reduce(mergewith(+), sf_list)

    @info "-------------- $tag shapes begin ------------ "
    shape_tasks = mapreduce(shape_variation -> WVZAnalysis.prep_tasks(tag; shape_variation), vcat,
                                 WVZAnalysis.SHAPE_TREE_NAMES)
    println("$(length(shape_tasks)) tasks in total")
    shape_list = @showprogress pmap(main_looper, shape_tasks)
    shape_hist = reduce(mergewith(+), shape_list)

    Hs = merge(sf_hist, shape_hist)
    serialize(joinpath(p, "$(tag)_pmap.jlserialize"), Hs)
    Hs
end

