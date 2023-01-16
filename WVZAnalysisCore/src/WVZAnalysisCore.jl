module WVZAnalysisCore

using UnROOT, FHist, LinearAlgebra, LorentzVectorHEP, LazyArrays, JSON3, ProgressMeter
using XGBoost, Distributed, Serialization
using ONNX, ONNX.Umlaut 

export AnalysisTask
export prep_tasks, main_looper, hist_root, arrow_making

include("./analysis_utils.jl")
include("./constants.jl")
include("./WWZ_ana.jl")
include("./arrow_exfiltration.jl")
include("./hist_exfiltration.jl")
include("./mainlooper.jl")

end
