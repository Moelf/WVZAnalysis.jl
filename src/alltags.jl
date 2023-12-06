const SIG_TAGS = ("Signal_WWZ", "Signal_WZZ", "Signal_ZZZ", "Signal_ZH_ZWW", "Signal_WH_WZZ", "Signal_ZH_ZZZ")
const BKG_TAGS = ("ZZ0j", "ZZ1j", "ZZ2plusj",
 # "Zjets",
 # "ttbar",
 # "WZ",
 # "tZ",
 "ttZ",
 "tWZ",
 "VBS",
 "WH_Bkg",
 "ZH_Bkg",
 "Others")
const ALL_TAGS = [SIG_TAGS...; BKG_TAGS...]
