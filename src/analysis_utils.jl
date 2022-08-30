"""
extract dsid from a file name, used to match with systematic files
"""
const MINITREE_DIR = Ref("/data/jiling/WVZ/v2.3")

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

    prog = Progress(mapreduce(length∘readdir, +, dirs), 0.5)
    isdata = (lowercase(tag) == "data")
    files = mapreduce(vcat, dirs) do d
        sfsys_dir(d; prog, isdata)
    end
    sort!(files; by=x->filesize(x[2]), rev=true)
    println("$tag starting:")
    ex = WorkStealingEx(; basesize=length(files) ÷ 4)
    @floop ex for (sumWeight, F) in files
        x = WVZAnalysis.main_looper(F; treename, sfsyst=false, sumWeight, isdata)
        if prog !== nothing 
            next!(prog)
        end
        @reduce(s .+= x)
    end
    return s
end

function sfsys(tag; controlregion = nothing)
    dirs = root_dirs(tag; variation = "sf")

    prog = Progress(mapreduce(length∘readdir, +, dirs), 0.5)
    isdata = (lowercase(tag) == "data")
    files = mapreduce(vcat, dirs) do d
        sfsys_dir(d; prog, isdata)
    end
    sort!(files; by=x->filesize(x[2]), rev=true)
    println("$tag starting:")
    # ex = ThreadedEx(; basesize = 1)
    ex = WorkStealingEx(; basesize=length(files) ÷ 4)
    @floop ex for (sumWeight, F) in files
        x = WVZAnalysis.main_looper(F; sfsyst=false, sumWeight, isdata, controlregion)
        if prog !== nothing 
            next!(prog)
        end
        @reduce(s .+= x)
    end
    return s
end

function sfsys_dir(dir_path; prog = nothing, isdata=false)
    files = filter!(endswith(".root"), readdir(dir_path; join = true))
    sumsum = sumsumWeight(files)
    return (sumsum .=> files)
end

function arrow_making(tag)
    dirs = root_dirs(tag; variation = "sf")
    isdata = lowercase(tag) == "data"
    prog = Progress(mapreduce(length∘readdir, +, dirs), 0.3)
    println("$tag starting:")
    res = ThreadsX.map(dirs) do d
        arrow_making_dir(d; prog, isdata)
    end
    foldl((x,y) -> append!.(x,y), res)
    first(res)
end

function arrow_making_dir(dir_path; prog = nothing, isdata)
    files = filter!(endswith(".root"), readdir(dir_path; join = true))
    sumWeight = sumsumWeight(files)
    res = map(files) do F
        r = WVZAnalysis.main_looper(F; sfsyst=false, sumWeight, arrow_making=true, isdata)
        if prog !== nothing 
            next!(prog)
        end
        r
    end
    foldl((x,y) -> append!.(x,y), res)
    first(res)
end
