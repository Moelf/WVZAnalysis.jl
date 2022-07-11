"""
extract dsid from a file name, used to match with systematic files
"""
const MINITREE_DIR = Ref("/data/jiling/WVZ/v2.2")

function extrac_dsid(str)
    return match(r"(\d{6})", str).captures[1] #dsid
    # match(r"e\d{4}_.*r\d+", str).match #amitag
end

"""
    root_dirs(tag::AbstractString; variation = "sf")

Reutrn a list of dir paths associated with a `tag` and variation.
"""
function root_dirs(tag::AbstractString; variation = "sf")
    _LIST = joinpath(dirname(@__DIR__), "lists/v2.2.list.json") |> read |> JSON3.read
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
            any(occursin(folder_name), dsids) && contains(folder_name, "v2.2"*variation)
        end
        return joinpath.(MINITREE_DIR[], sel2)
    end
end

"""
multiple tags
"""
root_dirs(tags; variation = "sf") = mapreduce(x->root_dirs(x; variation), vcat, unique(tags))

function sumsumWeight(paths)
    res = 0.0
    for p in paths
        r = ROOTFile(p)
        if !haskey(r, "sumWeight")
            return 1.0
        end
        res += r["sumWeight"][:fN][3]
    end
    return res
end

function shapesys(tag, treename)
    dirs = root_dirs(tag; variation = "shape")
    println("$tag-$treename starting:")
    prog = Progress(mapreduce(length∘readdir, +, dirs), 0.3)
    mapreduce((.+), dirs) do d
        shapesys_dir(d, treename; prog)
    end
end

function shapesys_dir(dir_path, treename; prog = nothing)
    files = filter!(endswith(".root"), readdir(dir_path; join = true))
    sumWeight = sumsumWeight(files)
    return mapreduce((.+), files) do F
        r = WVZAnalysis.main_looper(F; treename, sumWeight)
        if prog !== nothing 
            next!(prog)
        end
        r
    end
end

function sfsys(tag)
    dirs = root_dirs(tag; variation = "sf")
    prog = Progress(mapreduce(length∘readdir, +, dirs), 0.3)
    println("$tag starting:")
    isdata = (lowercase(tag) == "data")
    mapreduce((.+), dirs) do d
        sfsys_dir(d; prog, isdata)
    end
end

function sfsys_dir(dir_path; prog = nothing, isdata=false)
    files = filter!(endswith(".root"), readdir(dir_path; join = true))
    sumWeight = sumsumWeight(files)
    return mapreduce((.+), files) do F
        r = WVZAnalysis.main_looper(F; sfsyst=false, sumWeight, isdata)
        if prog !== nothing 
            next!(prog)
        end
        r
    end
end

function arrow_making(tag)
    dirs = root_dirs(tag; variation = "sf")
    prog = Progress(mapreduce(length∘readdir, +, dirs), 0.3)
    println("$tag starting:")
    res = ThreadsX.map(dirs) do d
        arrow_making_dir(d; prog)
    end
    foldl((x,y) -> append!.(x,y), res)
    first(res)
end

function arrow_making_dir(dir_path; prog = nothing)
    files = filter!(endswith(".root"), readdir(dir_path; join = true))
    sumWeight = sumsumWeight(files)
    res = map(files) do F
        r = WVZAnalysis.main_looper(F; sfsyst=false, sumWeight, arrow_making=true)
        if prog !== nothing 
            next!(prog)
        end
        r
    end
    foldl((x,y) -> append!.(x,y), res)
    first(res)
end
