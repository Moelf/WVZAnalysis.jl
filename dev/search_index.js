var documenterSearchIndex = {"docs":
[{"location":"internalapis/#Internal","page":"Internal APIs","title":"Internal","text":"","category":"section"},{"location":"internalapis/","page":"Internal APIs","title":"Internal APIs","text":"Modules = [WVZAnalysis, WVZReportExt]\nFilter = t -> t ∉ (WVZAnalysis.AnalysisTask,\n    WVZAnalysis.prep_tasks,\n    WVZAnalysis.main_looper,\n    WVZAnalysis.arrow_main)","category":"page"},{"location":"internalapis/#WVZAnalysis.arrow_making-Tuple{Any}","page":"Internal APIs","title":"WVZAnalysis.arrow_making","text":"arrow_making(tasks; mapper = map)\n\nTake a collection of tasks, run them via map and mergewith(append!). Returns a dict of vectors representing the datas after filtering.\n\n\n\n\n\n","category":"method"},{"location":"internalapis/#WVZAnalysis.hist_main-Tuple{Any}","page":"Internal APIs","title":"WVZAnalysis.hist_main","text":"hist_main(tag; mapfun=robust_pmap, no_shape = false, output_dir, kw...)\n\nRunning the main looper for a tag (e.g. VH, ttZ) to produce histogram kw... takes anything that AnalysisTask takes.\n\n\n\n\n\n","category":"method"},{"location":"internalapis/#WVZAnalysis.make_sfsys_wgt!","page":"Internal APIs","title":"WVZAnalysis.make_sfsys_wgt!","text":"make_sfsys_wgt!(evt, wgt, nominal_branch, wgt_name, wgt_idx=1)\n\nuse wgt_name to find alternative branch names\nFor every variation, make a new entry in wgt::Dict.\nAlso multiply nominal for wgt[:NOMINAL]\n\n\n\n\n\n","category":"function"},{"location":"internalapis/#WVZAnalysis.root_dirs-Tuple{AbstractString}","page":"Internal APIs","title":"WVZAnalysis.root_dirs","text":"root_dirs(tag::AbstractString; variation = \"sf\")\n\nReutrn a list of dir paths associated with a tag +sf or +shape\n\n\n\n\n\n","category":"method"},{"location":"internalapis/#WVZAnalysis.set_bdt_model_dir-Tuple{String}","page":"Internal APIs","title":"WVZAnalysis.set_bdt_model_dir","text":"Set where BDT modesl are (XGBoost)\n\n\n\n\n\n","category":"method"},{"location":"internalapis/#WVZAnalysis.set_minitree_dir-Tuple{String}","page":"Internal APIs","title":"WVZAnalysis.set_minitree_dir","text":"Set where minitree is.\n\n\n\n\n\n","category":"method"},{"location":"internalapis/#WVZAnalysis.sumsumWeight-Tuple{Any}","page":"Internal APIs","title":"WVZAnalysis.sumsumWeight","text":"Compute the sum of weights given path to a directory, and cache the result in a text file.\n\n\n\n\n\n","category":"method"},{"location":"internalapis/#WVZAnalysis.@fill_dict!-Tuple{Any, Any, Any}","page":"Internal APIs","title":"WVZAnalysis.@fill_dict!","text":"Example:\n\njulia> dd = Dict(:Nlep => [], :lep1_pid=>[]);\n\njulia> Nlep = 10;\n\njulia> lep1_pid =11;\n\njulia> @fill_dict! dd push! Nlep,lep1_pid;\n\n# this expands to\npush!(dd[:Nlep], Nlep)\npush!(dd[:lep1_pid], lep1_pid)\n\njulia> dd\nDict{Symbol, Vector{Any}} with 2 entries:\n  :Nlep     => [10]\n  :lep1_pid => [11]\n\n\n\n\n\n","category":"macro"},{"location":"internalapis/#WVZReportExt.print_sigtable-Tuple{Any}","page":"Internal APIs","title":"WVZReportExt.print_sigtable","text":"print_sigtable(full_table; io=stdout)\n\nTakes the output of significance_table and pretty print it:\n\nExample\n\njulia> M = significance_table(<path>);\n\njulia> print_sigtable(M)\n┌──────────────┬────────────────┬───────────────┬──────────────┬──────────────┬───────────────┐\n│              │     SF-inZ     │    SF-noZ     │      DF      │    CR-ZZ     │    CR-ttZ     │\n├──────────────┼────────────────┼───────────────┼──────────────┼──────────────┼───────────────┤\n│    Signal    │  10.66 ± 0.07  │  9.31 ± 0.1   │ 10.73 ± 0.14 │ 1.07 ± 0.02  │  2.15 ± 0.04  │\n│      ZZ      │ 1219.92 ± 5.38 │ 469.06 ± 2.44 │ 19.78 ± 0.45 │ 555.68 ± 2.4 │ 69.19 ± 0.74  │\n│    Zjets     │  -0.02 ± 0.13  │  2.6 ± 2.22   │ 6.47 ± 5.51  │  -0.0 ± 0.0  │  0.63 ± 0.39  │\n│    Zgamma    │   0.0 ± 0.0    │   0.0 ± 0.0   │  0.3 ± 0.29  │  0.0 ± 0.0   │   0.0 ± 0.0   │\n│    ttbar     │   0.0 ± 0.0    │  0.63 ± 0.18  │  0.28 ± 0.1  │  0.0 ± 0.0   │  0.43 ± 0.13  │\n│      WZ      │   0.36 ± 0.1   │  1.79 ± 0.23  │ 2.24 ± 0.29  │  0.0 ± 0.0   │  0.29 ± 0.06  │\n│      tZ      │  0.01 ± 0.01   │  0.07 ± 0.03  │ 0.06 ± 0.02  │  0.0 ± 0.0   │  0.18 ± 0.04  │\n│     ttZ      │  1.23 ± 0.08   │  4.71 ± 0.16  │ 5.74 ± 0.18  │ 0.01 ± 0.01  │ 73.62 ± 0.61  │\n│     tWZ      │  0.56 ± 0.11   │  2.16 ± 0.23  │  2.5 ± 0.24  │ 0.01 ± 0.01  │ 14.71 ± 0.56  │\n│     VBS      │  10.23 ± 0.09  │  6.4 ± 0.08   │ 0.18 ± 0.01  │ 1.48 ± 0.03  │  0.83 ± 0.02  │\n│      VH      │  1.29 ± 0.71   │  5.76 ± 1.4   │ 5.77 ± 1.29  │  0.0 ± 0.0   │   0.0 ± 0.0   │\n├──────────────┼────────────────┼───────────────┼──────────────┼──────────────┼───────────────┤\n│    Others    │  0.06 ± 0.01   │  0.4 ± 0.13   │ 0.56 ± 0.08  │  0.0 ± 0.0   │  5.02 ± 0.09  │\n├──────────────┼────────────────┼───────────────┼──────────────┼──────────────┼───────────────┤\n│   Bkg Tot.   │ 1233.64 ± 5.43 │ 493.58 ± 3.61 │ 43.89 ± 5.7  │ 557.19 ± 2.4 │ 164.91 ± 1.19 │\n│ Significance │   0.3 ± 0.0    │  0.42 ± 0.0   │  1.56 ± 0.1  │  0.05 ± 0.0  │  0.17 ± 0.0   │\n└──────────────┴────────────────┴───────────────┴──────────────┴──────────────┴───────────────┘\n\n\n\n\n\n","category":"method"},{"location":"internalapis/#WVZReportExt.rebinscan-Tuple{Any, Any}","page":"Internal APIs","title":"WVZReportExt.rebinscan","text":"rebinscan(S, B; atleast=1, from=:right, by = (s,b) -> s/sqrt(b))\n\nTake two Hist1D, try to find the best bin-edges for rebinning, given these two conditions:\n\na metric to maximize, s and b are signal counts and background counts in each new bin\natleast specifcy the minimal number of s+b in each newly formed bin\n\nfrom keyword argument can be used to specifcy the direction of the scan, defautls to :right assuming the high-er end of the histogram is signal-like and one of the histograms is monotonic.\n\nThe function returns the values of new bin-edges including both ends (of course, the ends are identical to before).\n\nThe metric can also be replaced by the more accurate sqrt(2*((s + b) * log(1 + s/b) - s )).\n\n\n\n\n\n","category":"method"},{"location":"internalapis/#WVZReportExt.significance_table-Tuple{AbstractString}","page":"Internal APIs","title":"WVZReportExt.significance_table","text":"significance_table()\nsignificance_table(significance_matrix(); )\n\nExamples\n\njulia> M = significance_table()\n14×6 Matrix{Any}:\n \"Signal\"         10.66±0.067   …    1.072±0.017     2.151±0.037\n \"ZZ\"            1219.9±5.4          555.7±2.4       69.19±0.74\n \"Zjets\"         -0.019±0.13       -0.0004±0.0004     0.63±0.39\n \"Zgamma\"           0.0±0.0            0.0±0.0         0.0±0.0\n \"ttbar\"            0.0±0.0            0.0±0.0        0.43±0.13\n \"WZ\"             0.358±0.096   …   0.0044±0.0032    0.295±0.061\n \"tZ\"             0.012±0.012          0.0±0.0       0.179±0.043\n \"ttZ\"            1.231±0.08        0.0111±0.0074    73.62±0.61\n \"tWZ\"             0.56±0.11        0.0095±0.0095    14.71±0.56\n \"VBS\"           10.231±0.085        1.478±0.033     0.833±0.024\n \"VH\"              1.29±0.71    …      0.0±0.0         0.0±0.0\n \"Others\"        0.0568±0.0088      0.0021±0.0017     5.02±0.088\n \"Bkg Tot.\"      1233.6±5.4          557.2±2.4       164.9±1.2\n \"Significance\"  0.3031±0.002       0.0454±0.00072  0.1672±0.0029\n\n\n\n\n\n","category":"method"},{"location":"walkthrough/#Introduction-and-pre-requisite","page":"Step-by-step Walkthrough","title":"Introduction and pre-requisite","text":"","category":"section"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"The entire WVZAnalysis (4 lepton channel) can be break into rought 5 parts:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Make flat .arrow files from \"nominal\" samples\nUse .arrow files to train XGBoost models\nUse XGBoost models to produce histograms that include all systematics variations (saved temporarily via Serialization. This step involves HTCondor (or maybe other job system)\nConvert histograms to .root files via Python uproot\nRun fit with TRex-fitter (outside of this walkthrough)","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"This Julia package is almost self-contained, certainly all of Julia parts should be exactly reproducible. But we will rely on LCG release for the three Python pkg (uproot, hist, numpy), and because Python importing rules, at step 4 you need to start Julia with an environment variable, as noted below.","category":"page"},{"location":"walkthrough/#Before-you-start","page":"Step-by-step Walkthrough","title":"Before you start","text":"","category":"section"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"We will assume you have access to LCG release via CVMFS, and use the dev4 branch which contains Julia 1.10+.","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"You should also clone this repository and checkout the refactor_Julia_v1p9 branch (or if this branch is merged, just whatever release tag you need):","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"git clone https://github.com/Moelf/WVZAnalysis.jl/\ncd WVZAnalysis.jl","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Finally, instantiate the exact Julia versions we used for this analysis:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"# this is generally the command we need to use to start a Julia REPL\n# `.` here is the root directory of this git repo\nJULIA_CPU_TARGET=\"generic;znver2,clone_all\" JULIA_CONDAPKG_BACKEND=Null julia --project=. \n\n]instantiate # when you press `]` the prompt should switch to `(WVZAnalysis) pkg>`","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"note: Note\nIf you want to use LCG for all the Python stuff, remember to add this JULIA_CONDAPKG_BACKEND=Null which tells PythonCall.jl to not manage our Python environment!","category":"page"},{"location":"walkthrough/#Step-1","page":"Step-by-step Walkthrough","title":"Step 1","text":"","category":"section"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"using WVZAnalysis\n\nforeach(ALL_TAGS) do tag\n    arrow_main(tag; output_dir=\"/data/jiling/WVZ/v2.4.2-2024_06_06_arrow/\");\nend","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Remember to change output_dir to somewhere you have write access.","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"note: Note\nThis assumes WVZAnalysis.MINITREE_DIR points to the data dir, on AF UChicago, all the minitree files are copied to /data/jiling/WVZ/v2.4.2, if you need to change this, run WVZAnalysis.set_minitree_dir and restart Julia session. This setting presists across Julia sessions via Preference.jl, the values are stored in the LocalPreferences.toml file.","category":"page"},{"location":"walkthrough/#Step-2","page":"Step-by-step Walkthrough","title":"Step 2","text":"","category":"section"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Now that the arrow files are in place, we can train XGBoost, you really want to use GPU for this training, in the walkthrough, we will use AF UChicago's gpu node.","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"First we start an interactive condor job:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"condor_submit -interactive config/interactive.sub","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Then once you're dropped onto remote worker, you need to re setup environment:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"source /cvmfs/sft.cern.ch/lcg/views/dev4/latest/x86_64-centos7-gcc11-opt/setup.sh\n\ncd <where you cloned WVZAnalysis.jl>\n\nLD_LIBRARY_PATH='' julia --project=./WVZXGBoostExt","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Then, we need to instantiate on this node because hardware has changed (we have a GPU now):","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"]instantiate","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Finally, we can run our XGBoost training (remember to replace the paths):","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"using WVZXGBoostExt\n\ndf_all = WVZXGBoostExt.load_all_arrow(\"/data/jiling/WVZ/v2.4.2-2024_06_06_arrow/\")\nWVZXGBoostExt.train_and_log(df_all; output_dir = \"/data/jiling/WVZ/v2.4.2.-2024_06_12_hists/\", tree_method=\"gpu_hist\")","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"note: Note\nThe +queue=\"short\" and request_gpus=1 (in interactive.sub) are used by AF UChicago.","category":"page"},{"location":"walkthrough/#Step-3","page":"Step-by-step Walkthrough","title":"Step 3","text":"","category":"section"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"In this step, we use the cluster (HTCondor in this case) to run through all systematic variations:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"But first, you need to point our package to the new location that stores the XGBoost models:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"$ JULIA_CPU_TARGET=\"generic;znver2,clone_all\" JULIA_CONDAPKG_BACKEND=Null julia --project=. -e 'using WVZAnalysis; WVZAnalysis.set_bdt_model_dir(\"/data/jiling/WVZ/v2.3-2023_06_15_hists/\")'","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Then start a new Julia session:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"$ JULIA_CPU_TARGET=\"generic;znver2,clone_all\" JULIA_CONDAPKG_BACKEND=Null julia --project=.","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Then, we can schedule all the histogram making jobs:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"using ClusterManagers, Distributed, WVZAnalysis\n\naddprocs(HTCManager(80); extrajdl=[\"+queue=\\\"short\\\"\", \"request_memory = 4GB\"], exeflags = `--project=$(Base.active_project())`);\n\n@everywhere using WVZAnalysis","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"note: Note\nWe now request 80 jobs, and we no longer need gpu requirement since we're not training.  The init.jl fixes an issue where remote Julia workers sometimes can't find the correct network interface.","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"foreach([ALL_TAGS; \"Data\"]) do tag\n    hist_main(tag; output_dir=\"/data/jiling/WVZ/v2.3-2023_06_15_hists/\");\nend","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"You want to exit this Julia session after you're finished, which will also kill the condor jobs (unless you want to manually run rmprocs())","category":"page"},{"location":"walkthrough/#Step-4","page":"Step-by-step Walkthrough","title":"Step 4","text":"","category":"section"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"In this last step, we will convert the output from last step into .root files, start a new Julia session if you exit the last one:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"$ JULIA_CPU_TARGET=\"generic;znver2,clone_all\" JULIA_CONDAPKG_BACKEND=Null julia --project=.","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"and simply use the little Python \"extention\" package we wrote:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"julia> using WVZPythonExt\n\njulia> serial_to_root(\"/data/jiling/WVZ/v2.3-2023_06_15_hists/\")\n\njulia> filter(endswith(\"root\"), readdir(\"/data/jiling/WVZ/v2.3-2023_06_12_hists/\"))\n13-element Vector{String}:\n \"Data.root\"\n \"Others.root\"\n \"Signal.root\"\n \"VBS.root\"\n \"VH.root\"\n \"WZ.root\"\n \"ZZ.root\"\n \"Zgamma.root\"\n \"Zjets.root\"\n \"tWZ.root\"\n \"tZ.root\"\n \"ttZ.root\"\n \"ttbar.root\"","category":"page"},{"location":"#The-Big-Picture","page":"Introduction","title":"The Big Picture","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"like any Julia package, the source files are under /src\nmeta data files such as tag -> DSID mapping are under /config\nthe main looper is called main_looper() and all the real actions are in there\nthe atomic unit of runnable task is called AnalysisTask","category":"page"},{"location":"#A-typical-workflow","page":"Introduction","title":"A typical workflow","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Say we're trying to produce BDT score histogram for process ZZ:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"call prep_tasks(\"ZZ\") (this returns a vector of AnalysisTasks)\nrun main_looper over each of the tasks, you can use threading, or distributed parallelism.\neach main_looper would would return a dictionary of histograms\nmerge all the result (reduce(mergewith(+), results_from_3))\nprofit","category":"page"},{"location":"#Public-ish-APIs","page":"Introduction","title":"Public-ish APIs","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"These are the types and functions all a user would need,  in order to do pretty much anything.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The packages that are good to be familiar with: ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"UnROOT.jl\nFHist.jl","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"These functions/types and the direct callees are more documented than the rest of the library:","category":"page"},{"location":"#Analysis-and-data-dumping","page":"Introduction","title":"Analysis and data dumping","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"WVZAnalysis.AnalysisTask\nWVZAnalysis.prep_tasks\nWVZAnalysis.main_looper(t::AnalysisTask)\nWVZAnalysis.arrow_main","category":"page"},{"location":"#WVZAnalysis.AnalysisTask","page":"Introduction","title":"WVZAnalysis.AnalysisTask","text":"struct AnalysisTask\n    path::String\n    sumWeight::Float64\n    isdata::Bool = false\n    BDT_hist::Bool = true\n    arrow_making::Bool = false\n    sfsys::Bool = false\n    shape_variation::String = \"NOMINAL\"\n    controlregion::Symbol = :none\n    require_VHSig::Union{Nothing, Bool} = nothing\n    sel_Njet::Union{Nothing, Int} = nothing\nend\n\nFully define a task to be run on an executor, by calling main_looper(task::AnalysisTask); see also main_looper.\n\nYou most likely don't need to construct it manually, see prep_tasks:\n\nExample\n\njulia> prep_tasks(\"Signal\") |> first\npath=\"/data/jiling/WVZ/v2.3/user.jiling.WVZ_v2.3sf.363507.e6379_s3126_r10201_p4434_ANALYSIS.root/user.jiling.29896106._000001.ANALYSIS.root\"\nsumWeight=13812.79638671875\nisdata=false\nBDT_hist=true\narrow_making=false\nsfsys=false\nshape_variation=\"NOMINAL\"\ncontrolregion=:none\n\njulia> prep_tasks(\"Signal\"; arrow_making=true) |> first\nERROR: can't do produce arrow and NN histograms at the same time\nStacktrace:\n...\n..\n\njulia> prep_tasks(\"Signal\"; arrow_making=true, BDT_hist=false) |> first\npath=\"/data/jiling/WVZ/v2.3/user.jiling.WVZ_v2.3sf.363507.e6379_s3126_r10201_p4434_ANALYSIS.root/user.jiling.29896106._000001.ANALYSIS.root\"\nsumWeight=13812.79638671875\nisdata=false\nBDT_hist=false\narrow_making=true\nsfsys=false\nshape_variation=\"NOMINAL\"\ncontrolregion=:none\n\n\n\n\n\n","category":"type"},{"location":"#WVZAnalysis.prep_tasks","page":"Introduction","title":"WVZAnalysis.prep_tasks","text":"prep_tasks(tag; shape_variation=\"NOMINAL\", scouting=false, kw...)\n\nFor construction a collection of AnalysisTasks given tag and job options, check out [@AnalysisTask].\n\nAll posible tags can be found in config/file_list.json, there's also a convinient variable WVZAnalysis.ALL_TAGS that keeps track of all processes in an exclusive mannar:\n\nExample\n\njulia> WVZAnalysis.ALL_TAGS\n(\"Signal\", \"ZZ\", \"Zjets\", \"Zgamma\", \"ttbar\", \"WZ\", \"tZ\", \"ttZ\", \"tWZ\", \"VBS\", \"VH\", \"Others\")\n\njulia> all_nominal_tasks = mapreduce(prep_tasks, vcat, WVZAnalysis.ALL_TAGS);\n\njulia> length(all_nominal_tasks)\n768\n\n\n\n\n\n","category":"function"},{"location":"#WVZAnalysis.main_looper-Tuple{AnalysisTask}","page":"Introduction","title":"WVZAnalysis.main_looper","text":"main_looper(task::AnalysisTask)\n\nThe main entry point for running the main looper for a given task, it's done in steps:\n\ndestruct all the options from the task (see AnalysisTask):\n\n    (; path, sumWeight, arrow_making, BDT_hist, isdata, \n     shape_variation, controlregion, sfsys) = task\n\ndetermine what output dict to prepare.\ncall the lower level main_looper (which has all arguments explicitly laied out)\n\nThis function also serves a function barrier for performance reason, because we have so many different behaviors in the same loop function.\n\n\n\n\n\n","category":"method"},{"location":"#WVZAnalysis.arrow_main","page":"Introduction","title":"WVZAnalysis.arrow_main","text":"arrow_main(tag; mapfun=ThreadsX.map, output_dir, kw...)\n\nRunning the main looper for a tag (e.g. VH, ttZ) to produce arrow kw... takes anything that AnalysisTask takes.\n\n\n\n\n\n","category":"function"},{"location":"#Result-reporting","page":"Introduction","title":"Result reporting","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"WVZReportExt.significance_table\nWVZReportExt.print_sigtable","category":"page"},{"location":"#WVZReportExt.significance_table","page":"Introduction","title":"WVZReportExt.significance_table","text":"significance_table()\nsignificance_table(significance_matrix(); )\n\nExamples\n\njulia> M = significance_table()\n14×6 Matrix{Any}:\n \"Signal\"         10.66±0.067   …    1.072±0.017     2.151±0.037\n \"ZZ\"            1219.9±5.4          555.7±2.4       69.19±0.74\n \"Zjets\"         -0.019±0.13       -0.0004±0.0004     0.63±0.39\n \"Zgamma\"           0.0±0.0            0.0±0.0         0.0±0.0\n \"ttbar\"            0.0±0.0            0.0±0.0        0.43±0.13\n \"WZ\"             0.358±0.096   …   0.0044±0.0032    0.295±0.061\n \"tZ\"             0.012±0.012          0.0±0.0       0.179±0.043\n \"ttZ\"            1.231±0.08        0.0111±0.0074    73.62±0.61\n \"tWZ\"             0.56±0.11        0.0095±0.0095    14.71±0.56\n \"VBS\"           10.231±0.085        1.478±0.033     0.833±0.024\n \"VH\"              1.29±0.71    …      0.0±0.0         0.0±0.0\n \"Others\"        0.0568±0.0088      0.0021±0.0017     5.02±0.088\n \"Bkg Tot.\"      1233.6±5.4          557.2±2.4       164.9±1.2\n \"Significance\"  0.3031±0.002       0.0454±0.00072  0.1672±0.0029\n\n\n\n\n\n","category":"function"},{"location":"#WVZReportExt.print_sigtable","page":"Introduction","title":"WVZReportExt.print_sigtable","text":"print_sigtable(full_table; io=stdout)\n\nTakes the output of significance_table and pretty print it:\n\nExample\n\njulia> M = significance_table(<path>);\n\njulia> print_sigtable(M)\n┌──────────────┬────────────────┬───────────────┬──────────────┬──────────────┬───────────────┐\n│              │     SF-inZ     │    SF-noZ     │      DF      │    CR-ZZ     │    CR-ttZ     │\n├──────────────┼────────────────┼───────────────┼──────────────┼──────────────┼───────────────┤\n│    Signal    │  10.66 ± 0.07  │  9.31 ± 0.1   │ 10.73 ± 0.14 │ 1.07 ± 0.02  │  2.15 ± 0.04  │\n│      ZZ      │ 1219.92 ± 5.38 │ 469.06 ± 2.44 │ 19.78 ± 0.45 │ 555.68 ± 2.4 │ 69.19 ± 0.74  │\n│    Zjets     │  -0.02 ± 0.13  │  2.6 ± 2.22   │ 6.47 ± 5.51  │  -0.0 ± 0.0  │  0.63 ± 0.39  │\n│    Zgamma    │   0.0 ± 0.0    │   0.0 ± 0.0   │  0.3 ± 0.29  │  0.0 ± 0.0   │   0.0 ± 0.0   │\n│    ttbar     │   0.0 ± 0.0    │  0.63 ± 0.18  │  0.28 ± 0.1  │  0.0 ± 0.0   │  0.43 ± 0.13  │\n│      WZ      │   0.36 ± 0.1   │  1.79 ± 0.23  │ 2.24 ± 0.29  │  0.0 ± 0.0   │  0.29 ± 0.06  │\n│      tZ      │  0.01 ± 0.01   │  0.07 ± 0.03  │ 0.06 ± 0.02  │  0.0 ± 0.0   │  0.18 ± 0.04  │\n│     ttZ      │  1.23 ± 0.08   │  4.71 ± 0.16  │ 5.74 ± 0.18  │ 0.01 ± 0.01  │ 73.62 ± 0.61  │\n│     tWZ      │  0.56 ± 0.11   │  2.16 ± 0.23  │  2.5 ± 0.24  │ 0.01 ± 0.01  │ 14.71 ± 0.56  │\n│     VBS      │  10.23 ± 0.09  │  6.4 ± 0.08   │ 0.18 ± 0.01  │ 1.48 ± 0.03  │  0.83 ± 0.02  │\n│      VH      │  1.29 ± 0.71   │  5.76 ± 1.4   │ 5.77 ± 1.29  │  0.0 ± 0.0   │   0.0 ± 0.0   │\n├──────────────┼────────────────┼───────────────┼──────────────┼──────────────┼───────────────┤\n│    Others    │  0.06 ± 0.01   │  0.4 ± 0.13   │ 0.56 ± 0.08  │  0.0 ± 0.0   │  5.02 ± 0.09  │\n├──────────────┼────────────────┼───────────────┼──────────────┼──────────────┼───────────────┤\n│   Bkg Tot.   │ 1233.64 ± 5.43 │ 493.58 ± 3.61 │ 43.89 ± 5.7  │ 557.19 ± 2.4 │ 164.91 ± 1.19 │\n│ Significance │   0.3 ± 0.0    │  0.42 ± 0.0   │  1.56 ± 0.1  │  0.05 ± 0.0  │  0.17 ± 0.0   │\n└──────────────┴────────────────┴───────────────┴──────────────┴──────────────┴───────────────┘\n\n\n\n\n\n","category":"function"}]
}
