module WVZAnalysis

using UnROOT, FHist, LinearAlgebra, LorentzVectors, LoopVectorization, StaticArrays

include("./utils.jl")
include("./ZZZ_ana.jl")
include("./WWZ_ana.jl")
include("./WZZ_ana.jl")
include("./ana.jl")


end
