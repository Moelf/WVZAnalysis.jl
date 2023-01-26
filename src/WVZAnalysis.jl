module WVZAnalysis

const ANALYSIS_DIR = Ref("/data/jiling/WVZ/v2.3_hists")
const SIG_TAGS = ("Signal", )
const BKG_TAGS = ("ZZ", "Zjets", "Zgamma", "ttbar", "WZ", "tZ", "ttZ", "tWZ", "VBS", "VH", "Others")
const ALL_TAGS = [SIG_TAGS...; BKG_TAGS...]

if contains(get(ENV, "PYTHONHOME", ""), "cern.ch/lcg/release")
    if get(ENV, "JULIA_CONDAPKG_BACKEND", "") != "Null"
        @warn "Detected LCG Environment"
        @info "re-launch Julia with `env JULIA_CONDAPKG_BACKEND=Null julia --project=.`"

        error("Using LCG but Julia not launched with correct environment")
        # use LCG release python
    end
end

using FHist
using PythonCall, Serialization

export serial_to_root, arrow_making, hist_root, significance_table, print_sigtable, ANALYSIS_DIR

include("./reporting_utils.jl")
include("./python_utils.jl")

end
