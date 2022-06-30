module WVZAnalysis

using UnROOT, FHist, LinearAlgebra, LorentzVectorHEP,  Dictionaries

include("./utils.jl")
include("./ZZZ_ana.jl")
include("./WWZ_ana.jl")
include("./WZZ_ana.jl")
include("./ana.jl")
include("./mainlooper.jl")


end
