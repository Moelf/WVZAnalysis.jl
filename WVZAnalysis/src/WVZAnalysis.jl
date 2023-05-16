module WVZAnalysis

using Preferences

using UnROOT, FHist, LinearAlgebra, LorentzVectorHEP, SentinelArrays, JSON3, ProgressMeter
using XGBoost

export AnalysisTask, ALL_TAGS
export prep_tasks, main_looper, hist_root, arrow_making
export set_minitree_dir, set_bdt_model_dir

"""
    Set where minitree is.
"""
function set_minitree_dir(path::String)
    # Set it in our runtime values, as well as saving it to disk
    @set_preferences!("MINITREE_DIR" => path)
    @info("New path for MINITREE_DIR set; restart your Julia session for this change to take effect!")
end

"""
    Set where BDT modesl are (XGBoost)
"""
function set_bdt_model_dir(path::String)
    # Set it in our runtime values, as well as saving it to disk
    @set_preferences!("BDT_MODEL_DIR" => path)
    @info("New path for BDT_MODEL_DIR set; restart your Julia session for this change to take effect!")
end

const MINITREE_DIR = @load_preference("MINITREE_DIR", "/data/jiling/WVZ/v2.3")
const BDT_MODEL_DIR = @load_preference("BDT_MODEL_DIR", joinpath(dirname(@__DIR__), "BDT_models"))

include("./analysis_utils.jl")
include("./constants.jl")
include("./mainlooper.jl")

end
