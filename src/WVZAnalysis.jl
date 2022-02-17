module WVZAnalysis

using UnROOT, FHist, LinearAlgebra, LorentzVectors,  Dictionaries, Mixers

include("./utils.jl")
include("./ZZZ_ana.jl")
include("./WWZ_ana.jl")
include("./WZZ_ana.jl")
include("./ana.jl")
include("./arrow_exfiltration.jl")
include("./hist_exfiltration.jl")
include("./mainlooper.jl")

end
