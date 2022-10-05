module WVZAnalysis

using UnROOT, FHist, LinearAlgebra, LorentzVectorHEP, LazyArrays, JSON3, ProgressMeter
using ThreadsX, FoldsThreads, XGBoost
# using ONNX, ONNX.Ghost 

using Distributed

export sfsys, significance_table, print_sigtable, shapesys, arrow_making

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

end
