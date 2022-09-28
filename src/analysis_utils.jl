
"""
extract dsid from a file name, used to match with systematic files
"""
const MINITREE_DIR = Ref("/data/jiling/WVZ/v2.3")

const ONNX_MODEL_PATH = Ref("/data/grabanal/NN/NN_08_23.onnx")

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
    bst = Booster(; model_file = "/data/jiling/WVZ/v2.3-beta2_arrow/xgb_2022-09-27.model")
    return ary->(predict(bst, permutedims(ary))[1])::Float32
end


function NN_calc(model, scales, minimums, NN_input)::Float32
    @inbounds @simd for i in eachindex(NN_input, scales, minimums)
        NN_input[i] = fma(NN_input[i], scales[i], minimums[i])
    end
    return Ghost.play!(model, NN_input)[1]
end


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

Reutrn a list of dir paths associated with a `tag` and variation.
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

function sumsumWeight(paths)
    wgt_factor = if occursin(r"346645|346646|346647", paths[1])
        2.745e-4
    else
        1.0
    end
    return sumsumWeight(ROOTFile.(paths))*wgt_factor
end

function sumsumWeight(RFiles::Vector{ROOTFile})::Float64
    res = 0.0
    for r in RFiles
        if !haskey(r, "sumWeight")
            return 1.0
        end
        res += r["sumWeight"][:fN][3]
    end
    return res
end

function _runwork(RFiles, prog; kw...)
    println("processing $(length(RFiles)) root files in total.")
    s = ThreadsX.map(RFiles) do (sumWeight, R)
        next!(prog)
        x = WVZAnalysis.main_looper(R, sumWeight; kw...)
        x
    end
    finish!(prog)
    return foldl((.+), s)
end

function shapesys(tag, shape_variations::Vector; scouting=false, kw...)
    println("$tag ($shape_variations)")
    dirs = root_dirs(tag; variation = "shape")

    isdata = (lowercase(tag) == "data")
    if scouting
        dirs = first(dirs, 2)
    end

    prog = Progress(mapreduce(length∘readdir, +, dirs), 0.1)
    sum_and_RFiles = mapreduce(vcat, dirs) do d
        process_dir(d; scouting)
    end

    return _runwork(sum_and_RFiles, prog; shape_variation, isdata, kw...)
end

function sfsys(tag; scouting=false, kw...)
    println("$tag (SF)")
    dirs = root_dirs(tag; variation = "sf")

    isdata = (lowercase(tag) == "data")
    if scouting
        dirs = first(dirs, 2)
    end
    prog = Progress(mapreduce(length∘readdir, +, dirs), 0.1)
    files = mapreduce(vcat, dirs) do d
        process_dir(d; scouting)
    end
    return _runwork(files, prog; isdata, kw...)
end

function process_dir(dir_path; scouting = false)
    files = filter!(endswith(".root"), readdir(dir_path; join = true))
    if scouting
        files = first(files, 2)
    end
    # keep file handles
    RFiles = ROOTFile.(files)
    sumsum = sumsumWeight(files)
    return (sumsum .=> RFiles)
end

function arrow_making(tag)
    dirs = root_dirs(tag; variation = "sf")
    isdata = lowercase(tag) == "data"
    prog = Progress(mapreduce(length∘readdir, +, dirs), 0.5)
    println("$tag starting:")
    res = ThreadsX.map(dirs) do d
        arrow_making_dir(d; prog, isdata)
    end
    foldl((x,y) -> append!.(x,y), res)
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
