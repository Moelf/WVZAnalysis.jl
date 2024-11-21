module WVZAnalysis

using Preferences

using UnROOT, FHist, LinearAlgebra, LorentzVectorHEP, SentinelArrays, JSON3, ProgressMeter, Arrow
using XGBoost, Parallelism, Serialization, ThreadsX, OhMyThreads

import Parallelism
import ProgressMeter
ProgressMeter.ncalls(::typeof(Parallelism.robust_pmap), ::Function, args...) = ProgressMeter.ncalls_map(args...)


export AnalysisTask, ALL_TAGS
export prep_tasks, main_looper, hist_main, arrow_main
export set_minitree_dir, set_bdt_model_dir
export serial_to_root

using WVZReportExt
export significance_table, print_sigtable

include("./alltags.jl")
include("./constants.jl")
include("./analysis_utils.jl")
include("./mainlooper.jl")

function serial_to_root(p)
    error("`using PythonCall` first to trigger Extention load")
end

end
