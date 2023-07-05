var documenterSearchIndex = {"docs":
[{"location":"internalapis/#Internal","page":"Internal APIs","title":"Internal","text":"","category":"section"},{"location":"internalapis/","page":"Internal APIs","title":"Internal APIs","text":"Modules = [WVZAnalysis, WVZAnalysisCore]","category":"page"},{"location":"walkthrough/#Introduction-and-pre-requisite","page":"Step-by-step Walkthrough","title":"Introduction and pre-requisite","text":"","category":"section"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"The entire WVZAnalysis (4 lepton channel) can be break into rought 5 parts:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Make flat .arrow files from \"nominal\" samples\nUse .arrow files to train XGBoost models\nUse XGBoost models to produce histograms that include all systematics variations (saved temporarily via Serialization. This step involves HTCondor (or maybe other job system)\nConvert histograms to .root files via Python uproot\nRun fit with TRex-fitter (outside of this walkthrough)","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"This Julia package is almost self-contained, certainly all of Julia parts should be exactly reproducible. But we will rely on LCG release for the three Python pkg (uproot, hist, numpy), and because Python importing rules, at step 4 you need to start Julia with an environment variable, as noted below.","category":"page"},{"location":"walkthrough/#Before-you-start","page":"Step-by-step Walkthrough","title":"Before you start","text":"","category":"section"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"We will assume you have access to LCG release via CVMFS, and use the dev4 branch which contains Julia 1.9+:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"source /cvmfs/sft.cern.ch/lcg/views/dev4/latest/x86_64-centos7-gcc11-opt/setup.sh","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"You should also clone this repository and checkout the refactor_Julia_v1p9 branch (or if this branch is merged, just whatever release tag you need):","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"git clone https://github.com/Moelf/WVZAnalysis.jl/\ncd WVZAnalysis.jl\ngit checkout refactor_Julia_v1p9","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Finally, instantiate the exact Julia versions we used for this analysis:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"# this is generally the command we need to use to start a Julia REPL\n# `.` here is the root directory of this git repo\nJULIA_CPU_TARGET=generic JULIA_CONDAPKG_BACKEND=Null julia --project=. \n\n]instantiate # when you press `]` the prompt should switch to `(WVZAnalysis) pkg>`","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"note: Note\nIf you want to use LCG for all the Python stuff, remember to add this JULIA_CONDAPKG_BACKEND=Null which tells PythonCall.jl to not manage our Python environment!","category":"page"},{"location":"walkthrough/#Step-1","page":"Step-by-step Walkthrough","title":"Step 1","text":"","category":"section"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"using WVZAnalysis\n\nforeach(ALL_TAGS) do tag\n    arrow_main(tag; output_dir=\"/data/jiling/WVZ/v2.3-2023_06_06_arrow/\");\nend","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Remember to change output_dir to somewhere you have write access.","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"note: Note\nThis assumes WVZAnalysis.MINITREE_DIR points to the data dir, on AF UChicago, all the minitree files are copied to /data/jiling/WVZ/v2.3, if you need to change this, run WVZAnalysis.set_minitree_dir and restart Julia session. This setting presists across Julia sessions via Preference.jl, the values are stored in the LocalPreferences.toml file.","category":"page"},{"location":"walkthrough/#Step-2","page":"Step-by-step Walkthrough","title":"Step 2","text":"","category":"section"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Now that the arrow files are in place, we can train XGBoost, you really want to use GPU for this training, in the walkthrough, we will use AF UChicago's gpu node.","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"First we start an interactive condor job:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"condor_submit -interactive config/interactive.sub","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Then once you're dropped onto remote worker, you need to re setup environment:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"source /cvmfs/sft.cern.ch/lcg/views/dev4/latest/x86_64-centos7-gcc11-opt/setup.sh\n\ncd <where you cloned WVZAnalysis.jl>\n\nJULIA_CPU_TARGET=generic JULIA_CONDAPKG_BACKEND=Null julia --project=./WVZXGBoostExt","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Then, we need to instantiate on this node because hardware has changed (we have a GPU now):","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"]instantiate","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Finally, we can run our XGBoost training (remember to replace the paths):","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"using WVZXGBoostExt\n\ndf_all = WVZXGBoostExt.load_all_arrow(\"/data/jiling/WVZ/v2.3-2023_06_06_arrow/\")\nWVZXGBoostExt.train_and_log(df_all; output_dir = \"/data/jiling/WVZ/v2.3-2023_06_15_hists/\", tree_method=\"gpu_hist\")","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"note: Note\nThe +queue=\"short\" and request_gpus=1 (in interactive.sub) are used by AF UChicago.","category":"page"},{"location":"walkthrough/#Step-3","page":"Step-by-step Walkthrough","title":"Step 3","text":"","category":"section"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"In this step, we use the cluster (HTCondor in this case) to run through all systematic variations:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"But first, you need to point our package to the new location that stores the XGBoost models:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"$ JULIA_CPU_TARGET=generic JULIA_CONDAPKG_BACKEND=Null julia --project=. -e 'using WVZAnalysis; WVZAnalysis.set_bdt_model_dir(\"/data/jiling/WVZ/v2.3-2023_06_15_hists/\")'","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Then start a new Julia session:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"$ JULIA_CPU_TARGET=generic JULIA_CONDAPKG_BACKEND=Null julia --project=.","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"Then, we can schedule all the histogram making jobs:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"using ClusterManagers, Distributed, WVZAnalysis\n\naddprocs(HTCManager(100); extrajdl=[\"+queue=\\\"short\\\"\"], extraenv=[\"export JULIA_CPU_TARGET=generic\"], exeflags = `--project=$(Base.active_project()) -e 'include(\"/data/jiling/WVZ/init.jl\")'`);\n\n@everywhere using WVZAnalysis","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"note: Note\nWe now request 100 jobs, and we no longer need gpu requirement since we're not training.  The init.jl fixes an issue where remote Julia workers sometimes can't find the correct network interface.","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"foreach([ALL_TAGS; \"Data\"]) do tag\n    hist_main(tag; output_dir=\"/data/jiling/WVZ/v2.3-2023_06_15_hists/\");\nend","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"You want to exit this Julia session after you're finished, which will also kill the condor jobs (unless you want to manually run rmprocs())","category":"page"},{"location":"walkthrough/#Step-4","page":"Step-by-step Walkthrough","title":"Step 4","text":"","category":"section"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"In this last step, we will convert the output from last step into .root files, start a new Julia session if you exit the last one:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"$ JULIA_CPU_TARGET=generic JULIA_CONDAPKG_BACKEND=Null julia --project=.","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"and simply use the little Python \"extention\" package we wrote:","category":"page"},{"location":"walkthrough/","page":"Step-by-step Walkthrough","title":"Step-by-step Walkthrough","text":"julia> using WVZPythonExt\n\njulia> serial_to_root(\"/data/jiling/WVZ/v2.3-2023_06_15_hists/\")\n\njulia> filter(endswith(\"root\"), readdir(\"/data/jiling/WVZ/v2.3-2023_06_12_hists/\"))\n13-element Vector{String}:\n \"Data.root\"\n \"Others.root\"\n \"Signal.root\"\n \"VBS.root\"\n \"VH.root\"\n \"WZ.root\"\n \"ZZ.root\"\n \"Zgamma.root\"\n \"Zjets.root\"\n \"tWZ.root\"\n \"tZ.root\"\n \"ttZ.root\"\n \"ttbar.root\"","category":"page"},{"location":"#The-Big-Picture","page":"Introduction","title":"The Big Picture","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"like any Julia package, the source files are under /src\nmeta data files such as tag -> DSID mapping are under /config\nthe main looper is called main_looper() and all the real actions are in there\nthe atomic unit of runnable task is called AnalysisTask","category":"page"},{"location":"#A-typical-workflow","page":"Introduction","title":"A typical workflow","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Say we're trying to produce BDT score histogram for process ZZ:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"call prep_tasks(\"ZZ\") (this returns a vector of AnalysisTasks)\nrun main_looper over each of the tasks, you can use threading, or distributed parallelism.\neach main_looper would would return a dictionary of histograms\nmerge all the result (reduce(mergewith(+), results_from_3))\nprofit","category":"page"},{"location":"#Public-ish-APIs","page":"Introduction","title":"Public-ish APIs","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"These are the types and functions all a user would need,  in order to do pretty much anything.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The packages that are good to be familiar with: ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"UnROOT.jl\nFHist.jl","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"These functions/types and the direct callees are more documented than the rest of the library:","category":"page"},{"location":"#Analysis-and-data-dumping","page":"Introduction","title":"Analysis and data dumping","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"AnalysisTask\nprep_tasks\nmain_looper(t::AnalysisTask)\narrow_making\nhist_root","category":"page"},{"location":"#WVZAnalysis.AnalysisTask","page":"Introduction","title":"WVZAnalysis.AnalysisTask","text":"struct AnalysisTask\n    path::String\n    sumWeight::Float64\n    isdata::Bool = false\n    BDT_hist::Bool = true\n    arrow_making::Bool = false\n    sfsys::Bool = false\n    shape_variation::String = \"NOMINAL\"\n    controlregion::Symbol = :none\nend\n\nFully define a task to be run on an executor, by calling main_looper(task::AnalysisTask); see also main_looper.\n\nYou most likely don't need to construct it manually, see prep_tasks:\n\nExample\n\njulia> prep_tasks(\"Signal\") |> first\npath=\"/data/jiling/WVZ/v2.3/user.jiling.WVZ_v2.3sf.363507.e6379_s3126_r10201_p4434_ANALYSIS.root/user.jiling.29896106._000001.ANALYSIS.root\"\nsumWeight=13812.79638671875\nisdata=false\nBDT_hist=true\narrow_making=false\nsfsys=false\nshape_variation=\"NOMINAL\"\ncontrolregion=:none\n\njulia> prep_tasks(\"Signal\"; arrow_making=true) |> first\nERROR: can't do produce arrow and NN histograms at the same time\nStacktrace:\n...\n..\n\njulia> prep_tasks(\"Signal\"; arrow_making=true, BDT_hist=false) |> first\npath=\"/data/jiling/WVZ/v2.3/user.jiling.WVZ_v2.3sf.363507.e6379_s3126_r10201_p4434_ANALYSIS.root/user.jiling.29896106._000001.ANALYSIS.root\"\nsumWeight=13812.79638671875\nisdata=false\nBDT_hist=false\narrow_making=true\nsfsys=false\nshape_variation=\"NOMINAL\"\ncontrolregion=:none\n\n\n\n\n\n","category":"type"},{"location":"#WVZAnalysis.prep_tasks","page":"Introduction","title":"WVZAnalysis.prep_tasks","text":"prep_tasks(tag; shape_variation=\"NOMINAL\", scouting=false, kw...)\n\nFor construction a collection of AnalysisTasks given tag and job options, check out [@AnalysisTask].\n\nAll posible tags can be found in config/file_list.json, there's also a convinient variable WVZAnalysis.ALL_TAGS that keeps track of all processes in an exclusive mannar:\n\nExample\n\njulia> WVZAnalysis.ALL_TAGS\n(\"Signal\", \"ZZ\", \"Zjets\", \"Zgamma\", \"ttbar\", \"WZ\", \"tZ\", \"ttZ\", \"tWZ\", \"VBS\", \"VH\", \"Others\")\n\njulia> all_nominal_tasks = mapreduce(prep_tasks, vcat, WVZAnalysis.ALL_TAGS);\n\njulia> length(all_nominal_tasks)\n768\n\n\n\n\n\n","category":"function"},{"location":"#WVZAnalysis.main_looper-Tuple{AnalysisTask}","page":"Introduction","title":"WVZAnalysis.main_looper","text":"main_looper(task::AnalysisTask)\n\nThe main entry point for running the main looper for a given task, it's done in steps:\n\ndestruct all the options from the task (see AnalysisTask):\n\n    (; path, sumWeight, arrow_making, BDT_hist, isdata, \n     shape_variation, controlregion, sfsys) = task\n\ndetermine what output dict to prepare.\ncall the lower level main_looper (which has all arguments explicitly laied out)\n\nThis function also serves a function barrier for performance reason, because we have so many different behaviors in the same loop function.\n\n\n\n\n\n","category":"method"},{"location":"#Result-reporting","page":"Introduction","title":"Result reporting","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"significance_table\nprint_sigtable","category":"page"}]
}
