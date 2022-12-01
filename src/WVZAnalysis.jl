module WVZAnalysis

using UnROOT, FHist, LinearAlgebra, LorentzVectorHEP, LazyArrays, JSON3, ProgressMeter
using XGBoost
using ONNX, ONNX.Umlaut 

using Distributed

export AnalysisTask
export prep_tasks, main_looper, significance_table, print_sigtable, arrow_making, hist_root, hist_root_pmap

include("./analysis_utils.jl")
include("./constants.jl")
include("./ZZZ_ana.jl")
include("./WWZ_ana.jl")
include("./WZZ_ana.jl")
include("./ana.jl")
include("./arrow_exfiltration.jl")
include("./hist_exfiltration.jl")
include("./mainlooper.jl")
include("./reporting_utils.jl")
# include("./python_utils.jl")

end
