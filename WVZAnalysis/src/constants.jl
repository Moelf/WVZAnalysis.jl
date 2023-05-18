const CUTFLOW_NAMES = [
    :in_minitree,
    :pass_trigger,
    :fourl_with_eta_mask,
    :has_Z_candidate,
    :pass_lepton_quality,
    :pass_lepton_iso,
    :SFOS_m_12GeVplus,
    :lepton_minimum_pts,
    :dilepton_dR_0p1plus,
    :no_Bjet,
    :MET_10GeVplus
    ]
const SIG_TAGS = ("Signal", )
const BKG_TAGS = ("ZZ", "Zjets", "Zgamma", "ttbar", "WZ", "tZ", "ttZ", "tWZ", "VBS", "VH", "Others")
const ALL_TAGS = [SIG_TAGS...; BKG_TAGS...]
const Z_m = 91.1876 # everything in GeV
const e_mass = 0.51099885 / 1000
const m_mass = 105.65837 / 1000
function makeud(v)
    res = Symbol[]
    for each in v
        push!(res, Symbol(each, :__1up))
        push!(res, Symbol(each, :__1down))
    end
    return res
end
const e_pliv = [:EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR] |> makeud
const m_pliv = [:MUON_EFF_ISO_SYS, :MUON_EFF_ISO_STAT] |> makeud
const e_wgt_list = [:EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR] |> makeud
const e_recowgt_list = [:EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR] |> makeud
const m_wgt_list = [:MUON_EFF_RECO_STAT_LOWPT, :MUON_EFF_RECO_SYS, :MUON_EFF_RECO_SYS_LOWPT, :MUON_EFF_RECO_STAT] |> makeud
const m_ttva_list = [:MUON_EFF_TTVA_STAT, :MUON_EFF_TTVA_SYS] |> makeud
const btag_wgt_list = [
                       :FT_EFF_Eigen_C_1,
                       :FT_EFF_Eigen_B_0,
                       :FT_EFF_Eigen_Light_0,
                       :FT_EFF_Eigen_B_1,
                       :FT_EFF_Eigen_Light_3,
                       :FT_EFF_extrapolation,
                       :FT_EFF_extrapolation_from_charm,
                       :FT_EFF_Eigen_C_2,
                       :FT_EFF_Eigen_C_0,
                       :FT_EFF_Eigen_Light_2,
                       :FT_EFF_Eigen_Light_1,
                      ] |> makeud

const weight_sf = [:JET_JvtEfficiency, :PRW_DATASF] |> makeud
const SF_BRANCH_DICT = Dict(
                            :weight => weight_sf,
                            :v_m_wgtIso_PflowLoose_VarRad => m_pliv,
                            :v_e_wgtIso_PLImprovedTight_Medium => e_pliv,
                            :v_e_wgtIso_Loose_VarRad_LooseBLayer => e_pliv,
                            :v_m_wgtIso_PLImprovedTight => m_pliv,
                            :v_e_wgtReco  => e_recowgt_list,
                            :v_e_wgtLoose  => e_wgt_list,
                            :v_e_wgtMedium => e_wgt_list,
                            :v_m_wgtLoose  => m_wgt_list,
                            :v_m_wgtMedium => m_wgt_list,
                            :v_j_wgt_btag77 => btag_wgt_list,
                            :v_m_wgtTTVA => m_ttva_list,
                            :w_sf_fjvt => [:JET_fJvtEfficiency__1down, :JET_fJvtEfficiency__1up]
                           )
const SHAPE_TREE_NAMES = [
                          "EG_RESOLUTION_ALL__1down",
                          "EG_RESOLUTION_ALL__1up",
                          "EG_SCALE_AF2__1down",
                          "EG_SCALE_AF2__1up",
                          "EG_SCALE_ALL__1down",
                          "EG_SCALE_ALL__1up",
                          "JET_EtaIntercalibration_Modelling__1down",
                          "JET_EtaIntercalibration_Modelling__1up",
                          "JET_Pileup_OffsetMu__1down",
                          "JET_Pileup_OffsetMu__1up",
                          "JET_Pileup_OffsetNPV__1down",
                          "JET_Pileup_OffsetNPV__1up",
                          "JET_Pileup_PtTerm__1down",
                          "JET_Pileup_PtTerm__1up",
                          "JET_BJES_Response__1down",
                          "JET_BJES_Response__1up",
                          "JET_EffectiveNP_1__1down",
                          "JET_EffectiveNP_1__1up",
                          "JET_EffectiveNP_2__1down",
                          "JET_EffectiveNP_2__1up",
                          "JET_EffectiveNP_3__1down",
                          "JET_EffectiveNP_3__1up",
                          "JET_EffectiveNP_4__1down",
                          "JET_EffectiveNP_4__1up",
                          "JET_EffectiveNP_5__1down",
                          "JET_EffectiveNP_5__1up",
                          "JET_EffectiveNP_6__1down",
                          "JET_EffectiveNP_6__1up",
                          "JET_EffectiveNP_7__1down",
                          "JET_EffectiveNP_7__1up",
                          "JET_EffectiveNP_8restTerm__1down",
                          "JET_EffectiveNP_8restTerm__1up",
                          "JET_EtaIntercalibration_NonClosure_2018data__1down",
                          "JET_EtaIntercalibration_NonClosure_2018data__1up",
                          "JET_EtaIntercalibration_NonClosure_highE__1down",
                          "JET_EtaIntercalibration_NonClosure_highE__1up",
                          "JET_EtaIntercalibration_NonClosure_negEta__1down",
                          "JET_EtaIntercalibration_NonClosure_negEta__1up",
                          "JET_EtaIntercalibration_NonClosure_posEta__1down",
                          "JET_EtaIntercalibration_NonClosure_posEta__1up",
                          "JET_EtaIntercalibration_TotalStat__1down",
                          "JET_EtaIntercalibration_TotalStat__1up",
                          "JET_Flavor_Composition__1down",
                          "JET_Flavor_Composition__1up",
                          "JET_Flavor_Response__1down",
                          "JET_Flavor_Response__1up",
                          "JET_JER_DataVsMC_MC16__1down",
                          "JET_JER_DataVsMC_MC16__1up",
                          "JET_JER_EffectiveNP_1__1down",
                          "JET_JER_EffectiveNP_1__1up",
                          "JET_JER_EffectiveNP_10__1down",
                          "JET_JER_EffectiveNP_10__1up",
                          "JET_JER_EffectiveNP_11__1down",
                          "JET_JER_EffectiveNP_11__1up",
                          "JET_JER_EffectiveNP_12restTerm__1down",
                          "JET_JER_EffectiveNP_12restTerm__1up",
                          "JET_JER_EffectiveNP_2__1down",
                          "JET_JER_EffectiveNP_2__1up",
                          "JET_JER_EffectiveNP_3__1down",
                          "JET_JER_EffectiveNP_3__1up",
                          "JET_JER_EffectiveNP_4__1down",
                          "JET_JER_EffectiveNP_4__1up",
                          "JET_JER_EffectiveNP_5__1down",
                          "JET_JER_EffectiveNP_5__1up",
                          "JET_JER_EffectiveNP_6__1down",
                          "JET_JER_EffectiveNP_6__1up",
                          "JET_JER_EffectiveNP_7__1down",
                          "JET_JER_EffectiveNP_7__1up",
                          "JET_JER_EffectiveNP_8__1down",
                          "JET_JER_EffectiveNP_8__1up",
                          "JET_JER_EffectiveNP_9__1down",
                          "JET_JER_EffectiveNP_9__1up",
                          "JET_Pileup_RhoTopology__1down",
                          "JET_Pileup_RhoTopology__1up",
                          "JET_PunchThrough_MC16__1down",
                          "JET_PunchThrough_MC16__1up",
                          "JET_SingleParticle_HighPt__1down",
                          "JET_SingleParticle_HighPt__1up",
                          "MET_SoftTrk_ResoPara",
                          "MET_SoftTrk_ResoPerp",
                          "MET_SoftTrk_Scale__1down",
                          "MET_SoftTrk_Scale__1up",
                          "MUON_CB__1down",
                          "MUON_CB__1up",
                          "MUON_ID__1down",
                          "MUON_ID__1up",
                          "MUON_MS__1down",
                          "MUON_MS__1up",
                          "MUON_SAGITTA_DATASTAT__1down",
                          "MUON_SAGITTA_DATASTAT__1up",
                          "MUON_SAGITTA_RESBIAS__1down",
                          "MUON_SAGITTA_RESBIAS__1up",
                          "MUON_SCALE__1down",
                          "MUON_SCALE__1up",
                         ]

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

