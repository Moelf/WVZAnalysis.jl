# Introduction and pre-requisite

The entire WVZAnalysis (4 lepton channel) can be break into rought 5 parts:
1. Make flat `.arrow` files from "nominal" samples
2. Use `.arrow` files to train XGBoost models
3. Use XGBoost models to produce histograms that include all systematics variations (saved
   temporarily via `Serialization`. This step involves HTCondor (or maybe other job system)
4. Convert histograms to `.root` files via Python `uproot`
5. Run fit with TRex-fitter (outside of this walkthrough)

This Julia package is almost self-contained, certainly all of Julia parts should be exactly
reproducible. But we will rely on `LCG` release for the three Python pkg (`uproot, hist, numpy`),
and because Python importing rules, at step 4 you need to start Julia with an environment variable,
as noted below.


## Before you start
We will assume you have access to LCG release via CVMFS, and use the `dev4` branch which contains
Julia 1.9+:

```bash
source /cvmfs/sft.cern.ch/lcg/views/dev4/latest/x86_64-centos7-gcc11-opt/setup.sh
```

You should also clone this repository and checkout the `refactor_Julia_v1p9` branch (or if this
branch is merged, just whatever release tag you need):
```bash
git clone https://github.com/Moelf/WVZAnalysis.jl/
cd WVZAnalysis.jl
git checkout refactor_Julia_v1p9
```

Finally, instantiate the exact Julia versions we used for this analysis:
```
JULIA_CONDAPKG_BACKEND=Null julia --project=. # `.` here is the root directory of this git repo
]instantiate # when you press `]` the prompt should switch to `pkg>`
```

!!! note
    If you want to use LCG for all the Python stuff, remember to add this
    `JULIA_CONDAPKG_BACKEND=Null` which tells `PythonCall.jl` to not manage our Python environment!


# Step 1

```julia
using WVZAnalysis

foreach(ALL_TAGS) do tag
    arrow_main(tag; output_dir="/data/jiling/WVZ/v2.3-2023_06_06_arrow/");
end
```

Remember to change `output_dir` to somewhere you have write access.

!!! note
    This assumes `WVZAnalysis.MINITREE_DIR` points to the data dir, on AF UChicago, all the minitree
    files are copied to `/data/jiling/WVZ/v2.3`, if you need to change this, run
    `WVZAnalysis.set_minitree_dir` and restart Julia session. This setting presists across Julia
    sessions via `Preference.jl`, the values are stored in the `LocalPreferences.toml` file.

# Step 2
Now that the arrow files are in place, we can train XGBoost, you should start a new Julia process
for this (and exit after done, since this reserve a GPU node):
```julia
using ClusterManagers, Distributed

addprocs(HTCManager(1); extrajdl=["+queue=\"short\"", "request_gpus=1"], extraenv=["export JULIA_CPU_TARGET=generic"], exeflags = `--project=$(Base.active_project()) -e 'include("/data/jiling/WVZ/init.jl")'`);

@fetchfrom 2 begin
    @eval using WVZXGBoostExt
    df_all = WVZXGBoostExt.load_all_arrow("/data/jiling/WVZ/v2.3-2023_06_06_arrow/")
    WVZXGBoostExt.train_and_log(df_all; output_dir = "/data/jiling/WVZ/v2.3-2023_06_15_hists/", tree_method="gpu_hist")
end
```

As before, remember to reaplce `output_dir` and input path with what you have

!!! note
    The `+queue="short"` and `request_gpus=1` are specific to AF UChicago, and the inclusion of
    `init.jl` fixes an issue where remote Julia workers sometimes can't find the correct network
    interface.


# Step 3
In this step, we use the cluster (HTCondor in this case) to run through all systematic variations:

But first, you need to point our package to the new location that stores the XGBoost models:
```bash
$ JULIA_CONDAPKG_BACKEND=Null julia --project=. -e 'using WVZAnalysis; WVZAnalysis.set_bdt_model_dir("/data/jiling/WVZ/v2.3-2023_06_15_hists/")'
```

Then restart Julia session as instructed,

```julia
using ClusterManagers, Distributed, WVZAnalysis

addprocs(HTCManager(100); extrajdl=["+queue=\"short\""], extraenv=["export JULIA_CPU_TARGET=generic"], exeflags = `--project=$(Base.active_project()) -e 'include("/data/jiling/WVZ/init.jl")'`);

@everywhere using WVZAnalysis
```

!!! note
    The difference between this `addprocs` from the last one is that we now request 100 jobs, and we
    dropped the gpu requirement since we're not training.

```julia
foreach([ALL_TAGS; "Data"]) do tag
    hist_main(tag; output_dir="/data/jiling/WVZ/v2.3-2023_06_15_hists/");
end
```

You want to exit this Julia session after you're finished, which will also kill the condor jobs
(unless you want to manually run `rmprocs()`)

# Step 4

Finally, we will convert the output from last step into `.root` files, start a new Julia session if
you exit the last one:
```bash
JULIA_CONDAPKG_BACKEND=Null julia --project=.
```

and simply use the little Python "extention" package we wrote:
```julia
julia> using WVZPythonExt

julia> serial_to_root("/data/jiling/WVZ/v2.3-2023_06_15_hists/")

julia> filter(endswith("root"), readdir("/data/jiling/WVZ/v2.3-2023_06_12_hists/"))
13-element Vector{String}:
 "Data.root"
 "Others.root"
 "Signal.root"
 "VBS.root"
 "VH.root"
 "WZ.root"
 "ZZ.root"
 "Zgamma.root"
 "Zjets.root"
 "tWZ.root"
 "tZ.root"
 "ttZ.root"
 "ttbar.root"
```

