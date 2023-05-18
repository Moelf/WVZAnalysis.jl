# WVZAnalysis

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://moelf.github.io/WVZAnalysis.jl/dev)
[![DOI](https://zenodo.org/badge/388537649.svg)](https://zenodo.org/badge/latestdoi/388537649)


## Words for Devs
If you're familiar with basic Julia workflow, check out the [documentation](https://moelf.github.io/WVZAnalysis.jl/dev)
directly.

## New to Julia:
- Download Julia, we require 1.9 or above (for Pkg extension functionality)
    - if you have access to `/cvmfs`, you can also access Julia via LCG release:
```
source /cvmfs/sft.cern.ch/lcg/views/dev4/latest/x86_64-centos7-gcc11-opt/setup.sh
```

- `$ cd <path to this repo folder>`

- Launch Julia via `julia --project=.`, this is because
you want to use `Manifest.toml` to keep exact versions of packages including all dependencies.

- For the first time only, you should run `]instantiate` to install all the depeendencies.

## Folder structre:
- The source code of analysis algorithm and cutflows etc. are inside `WVZAnalysis/src`


## Example
```bash
$ pwd
/home/jiling/.julia/dev/WVZAnalysis/ 

$ julia --project=.
```
then in Julia REPL:
```julia
julia> using WVZAnalysis

julia> tasks = prep_tasks("VH");

julia> tasks[1]
path="/data/jiling/WVZ/v2.3/user.jiling.WVZ_v2.3sf.342284.e4246_s3126_r10201_p4434_ANALYSIS.root/user.jiling.29896151._000001.ANALYSIS.root"
sumWeight=100000.0
isdata=false
BDT_hist=true
arrow_making=false
sfsys=false
shape_variation="NOMINAL"
controlregion=:none


julia> main_looper(tasks[1])
Dict{Symbol, FHist.Hist1D} with 19 entries:
  :SFnoZ__BDT__NOMINAL  => Hist1D{Float64}, edges=0.0:0.01:1.0, integral=0.0
  :CutFlow              => Hist1D{Int64}, edges=1:20, integral=149
  :SFinZCutFlow         => Hist1D{Int64}, edges=1:20, integral=0
  :DF__BDT__NOMINAL     => Hist1D{Float64}, edges=0.0:0.01:1.0, integral=0.4336316428757472
  :DFCutFlowWgt         => Hist1D{Float64}, edges=1:20, integral=3.582834263543141
  :SFnoZCutFlow         => Hist1D{Int64}, edges=1:20, integral=2
  :DF__Njet__NOMINAL    => Hist1D{Float64}, edges=0:5, integral=0.4336316428757472
  :SFnoZ__MET__NOMINAL  => Hist1D{Float64}, edges=0:5:300, integral=0.0
  :SFinZ__MET__NOMINAL  => Hist1D{Float64}, edges=0:5:300, integral=0.0
  :CutFlowWgt           => Hist1D{Float64}, edges=1:20, integral=77.51041686776865
  :ZZCR__Njet__NOMINAL  => Hist1D{Float64}, edges=0:5, integral=0.0
  :SFnoZCutFlowWgt      => Hist1D{Float64}, edges=1:20, integral=0.9438591695752855
  :SFnoZ__Njet__NOMINAL => Hist1D{Float64}, edges=0:5, integral=0.0
  :DFCutFlow            => Hist1D{Int64}, edges=1:20, integral=8
  :SFinZCutFlowWgt      => Hist1D{Float64}, edges=1:20, integral=0.0
  :SFinZ__BDT__NOMINAL  => Hist1D{Float64}, edges=0.0:0.01:1.0, integral=0.0
  :DF__MET__NOMINAL     => Hist1D{Float64}, edges=0:5:300, integral=0.4336316428757472
  :ttZCR__Njet__NOMINAL => Hist1D{Float64}, edges=0:5, integral=0.0
  :SFinZ__Njet__NOMINAL => Hist1D{Float64}, edges=0:5, integral=0.0
```
