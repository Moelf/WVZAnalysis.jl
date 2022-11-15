# WVZAnalysis

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://moelf.github.io/WVZAnalysis.jl/dev)
[![DOI](https://zenodo.org/badge/388537649.svg)](https://zenodo.org/badge/latestdoi/388537649)


## Words for Devs
If you're familiar with basic Julia workflow, check out the [documentation](https://moelf.github.io/WVZAnalysis.jl/dev)
directly.

## New to Julia:
- Download Julia 1.7 or above, we recommend 1.8 or above
    - if you have access to `/cvmfs`, you can also access Julia via LCG release:
```
source /cvmfs/sft.cern.ch/lcg/views/dev4/latest/x86_64-centos7-gcc11-opt/setup.sh
```

- We recommend lauching Julia via `julia --project=<path to this repo folder>`, this is because
you may want to use `Manifest.toml` to keep exact versions of packages including all dependencies.

- You should `]instantiate` or `]up` if you manually pulled this repo or the first time colneing it.


## Real-time out-of-core processing
```julia
julia> using WVZAnalysis

julia> using ClusterManagers, Distributed

julia> addprocs(HTCManager(10); extrajdl=["+queue=\"short\""], exeflags = `-e 'include("/data/jiling/WVZ/init.jl")'`);
Waiting for 10 workers: 1 2 3 4 5 6 7 8 9 10 .

julia> foreach(WVZAnalysis.ALL_TAGS) do tag
           WVZAnalysis.hist_root_pmap(tag; scouting=false, output_dir="/data/jiling/WVZ/v2.3_hists_uproot_nov14_andMET/");
       end
```
