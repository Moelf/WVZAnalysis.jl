# WVZAnalysis

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://moelf.github.io/WVZAnalysis.jl/dev)
[![DOI](https://zenodo.org/badge/388537649.svg)](https://zenodo.org/badge/latestdoi/388537649)


## Words for Devs
If you're familiar with basic Julia workflow, check out the [documentation](https://moelf.github.io/WVZAnalysis.jl/dev)
directly.

## New to Julia:
- Download Julia, we recommend 1.8 or above
    - if you have access to `/cvmfs`, you can also access Julia via LCG release:
```
source /cvmfs/sft.cern.ch/lcg/views/dev4/latest/x86_64-centos7-gcc11-opt/setup.sh
```

- Launch Julia via `julia --project=<path to this repo folder>`, this is because
you want to use `Manifest.toml` to keep exact versions of packages including all dependencies.

- You should `]instantiate` or `]up` if you manually pulled this repo or the first time colneing it.

## Folder structre:
- The analysis algorithm and cutflows are implemengted inside `WVZAnalysisCore` package
- The top-level (i.e. `WVZAnalysis`) is used for reporting, plotting, PythonCall -> `.root` and
most importantly to provide a `Manifest.toml` environment for exact reproduction.
- `WVZAnalysisCore/config` contains trained XGBoost model files and JSON files that map process name
(i.e. "Signal", "ttZ") to a list of files (DSIDs).


## Example workflow
```bash
$ pwd
/home/jiling/.julia/dev/WVZAnalysis/ 

$ julia --project=.
```
then in Julia REPL:
```julia
julia> using WVZAnalysisCore, ClusterManagers, Distributed

julia> addprocs(HTCManager(10); extrajdl=["+queue=\"short\""], exeflags = `--project=$(Base.active_project()) -e 'include("/data/jiling/WVZ/init.jl")'`);
Waiting for 10 workers: 1 2 3 4 5 6 7 8 9 10 .

julia> @everywhere using WVZAnalysisCore

julia> foreach(WVZAnalysisCore.ALL_TAGS) do tag
           hist_root(tag; scouting=false, output_dir="/data/jiling/WVZ/v2.3_hists_uproot_jan9/");
       end
```
