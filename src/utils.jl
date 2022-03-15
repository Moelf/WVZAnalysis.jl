const Z_m = 91.1876 * 10^3 # in MeV

pt(lv::LorentzVectorCyl) = lv.pt
eta(lv::LorentzVectorCyl) = lv.eta
phi(lv::LorentzVectorCyl) = lv.phi
mass(lv::LorentzVectorCyl) = lv.mass
Base.zero(lv::T) where T<:LorentzVectorCyl = T(0,0,0,0)

Base.zero(lv::T) where T<:LorentzVector = T(0,0,0,0)
mass(lv::LorentzVector) = sqrt(dot(lv, lv))
pt(lv::LorentzVector) = sqrt(lv.x^2 + lv.y^2)
mt2(lv::LorentzVector) = lv.t^2 - lv.z^2
mt(lv::LorentzVector) = mt2(lv)<0 ? -sqrt(-mt2(lv)) : sqrt(mt2(lv))
mag(lv::LorentzVector) = sqrt(lv.x^2 + lv.y^2 + lv.z^2)
@inline function CosTheta(lv::LorentzVector)
    fZ = lv.z
    ptot = mag(lv)
    return ifelse(ptot == 0.0, 1.0, fZ / ptot)
end
function eta(lv::LorentzVector)
    cosTheta = CosTheta(lv)
    (cosTheta^2 < 1.0) && return -0.5 * log((1.0 - cosTheta) / (1.0 + cosTheta))
    fZ = lv.z
    iszero(fZ) && return 0.0
    # Warning("PseudoRapidity","transvers momentum = 0! return +/- 10e10");
    fZ > 0.0 && return 10e10
    return -10e10
end
function phi(lv::LorentzVector)
    return (lv.x == 0.0 && lv.y == 0.0) ? 0.0 : atan(lv.y, lv.x)
end

function phi_mpi_pi(x)
    twopi = 2pi
    while (x >= pi)
        x -= twopi
    end
    while (x < -pi)
        x += twopi
    end
    return x
end

function deltaR(lv1, lv2)
    deta = eta(lv1) - eta(lv2)
    dphi = phi_mpi_pi(phi(lv1) - phi(lv2))
    return sqrt(deta * deta + dphi * dphi)
end

const _EISOS = (
    :HighPtCaloOnly,
    :TightTrackOnly_VarRad,
    :TightTrackOnly_FixedRad,
    :Tight_VarRad,
    :Loose_VarRad,
)
const _MISOS = (
    :PflowTight_VarRad,
    :PflowTight_FixedRad,
    :PflowLoose_VarRad,
    :PflowLoose_FixedRad,
    :HighPtTrackOnly,
    :TightTrackOnly_VarRad,
    :TightTrackOnly_FixedRad,
    :Tight_VarRad,
    :Tight_FixedRad,
    :Loose_VarRad,
    :Loose_FixedRad,
)

#FIXME thanks HEP, replace this with a @generated
function get_Isos(evt)
    v_l_passIso = Vector{Bool}[]

    v_e_passIso_TightTrackOnly_VarRad = evt.v_e_passIso_TightTrackOnly_VarRad
    v_e_passIso_Tight_VarRad = evt.v_e_passIso_Tight_VarRad
    v_e_passIso_Loose_VarRad = evt.v_e_passIso_Loose_VarRad
    v_m_passIso_PflowTight_VarRad = evt.v_m_passIso_PflowTight_VarRad
    v_m_passIso_PflowLoose_VarRad = evt.v_m_passIso_PflowLoose_VarRad

    @inbounds for i in eachindex(v_e_passIso_Loose_VarRad)
        push!(
            v_l_passIso,
            Bool[
                v_e_passIso_TightTrackOnly_VarRad[i],
                v_e_passIso_Tight_VarRad[i],
                v_e_passIso_Loose_VarRad[i],
            ],
        )
    end
    @inbounds for i in eachindex(v_m_passIso_PflowTight_VarRad)
        push!(
            v_l_passIso,
            Bool[
                 v_m_passIso_PflowTight_VarRad[i],
                 v_m_passIso_PflowLoose_VarRad[i],
            ],
        )
    end

    return v_l_passIso
end

const e_wgt_list = [:EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR]
const m_wgt_list = [:MUON_EFF_RECO_STAT_LOWPT, :MUON_EFF_RECO_SYS, :MUON_EFF_RECO_SYS_LOWPT, :MUON_EFF_RECO_STAT]
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
                      ]

const SF_BRANCH_DICT = Dict(
     :v_e_wgtLoose  => e_wgt_list,
     :v_e_wgtMedium => e_wgt_list,
     :v_m_wgtLoose  => m_wgt_list,
     :v_m_wgtMedium => m_wgt_list,
     :v_j_wgt_btag77 => btag_wgt_list,
    )

const SHAPE_TREE_NAMES = [
                          "tree_EG_RESOLUTION_ALL__1down",
                          "tree_EG_RESOLUTION_ALL__1up",
                          "tree_EG_SCALE_AF2__1down",
                          "tree_EG_SCALE_AF2__1up",
                          "tree_EG_SCALE_ALL__1down",
                          "tree_EG_SCALE_ALL__1up",
                          "tree_JET_BJES_Response__1down",
                          "tree_JET_BJES_Response__1up",
                          "tree_JET_EffectiveNP_1__1down",
                          "tree_JET_EffectiveNP_1__1up",
                          "tree_JET_EffectiveNP_2__1down",
                          "tree_JET_EffectiveNP_2__1up",
                          "tree_JET_EffectiveNP_3__1down",
                          "tree_JET_EffectiveNP_3__1up",
                          "tree_JET_EffectiveNP_4__1down",
                          "tree_JET_EffectiveNP_4__1up",
                          "tree_JET_EffectiveNP_5__1down",
                          "tree_JET_EffectiveNP_5__1up",
                          "tree_JET_EffectiveNP_6__1down",
                          "tree_JET_EffectiveNP_6__1up",
                          "tree_JET_EffectiveNP_7__1down",
                          "tree_JET_EffectiveNP_7__1up",
                          "tree_JET_EffectiveNP_8restTerm__1down",
                          "tree_JET_EffectiveNP_8restTerm__1up",
                          "tree_JET_EtaIntercalibration_Modelling__1down",
                          "tree_JET_EtaIntercalibration_Modelling__1up",
                          "tree_JET_EtaIntercalibration_NonClosure_2018data__1down",
                          "tree_JET_EtaIntercalibration_NonClosure_2018data__1up",
                          "tree_JET_EtaIntercalibration_NonClosure_highE__1down",
                          "tree_JET_EtaIntercalibration_NonClosure_highE__1up",
                          "tree_JET_EtaIntercalibration_NonClosure_negEta__1down",
                          "tree_JET_EtaIntercalibration_NonClosure_negEta__1up",
                          "tree_JET_EtaIntercalibration_NonClosure_posEta__1down",
                          "tree_JET_EtaIntercalibration_NonClosure_posEta__1up",
                          "tree_JET_EtaIntercalibration_TotalStat__1down",
                          "tree_JET_EtaIntercalibration_TotalStat__1up",
                          "tree_JET_Flavor_Composition__1down",
                          "tree_JET_Flavor_Composition__1up",
                          "tree_JET_Flavor_Response__1down",
                          "tree_JET_Flavor_Response__1up",
                          "tree_JET_JER_DataVsMC_MC16__1down",
                          "tree_JET_JER_DataVsMC_MC16__1up",
                          "tree_JET_JER_EffectiveNP_1__1down",
                          "tree_JET_JER_EffectiveNP_1__1up",
                          "tree_JET_JER_EffectiveNP_10__1down",
                          "tree_JET_JER_EffectiveNP_10__1up",
                          "tree_JET_JER_EffectiveNP_11__1down",
                          "tree_JET_JER_EffectiveNP_11__1up",
                          "tree_JET_JER_EffectiveNP_12restTerm__1down",
                          "tree_JET_JER_EffectiveNP_12restTerm__1up",
                          "tree_JET_JER_EffectiveNP_2__1down",
                          "tree_JET_JER_EffectiveNP_2__1up",
                          "tree_JET_JER_EffectiveNP_3__1down",
                          "tree_JET_JER_EffectiveNP_3__1up",
                          "tree_JET_JER_EffectiveNP_4__1down",
                          "tree_JET_JER_EffectiveNP_4__1up",
                          "tree_JET_JER_EffectiveNP_5__1down",
                          "tree_JET_JER_EffectiveNP_5__1up",
                          "tree_JET_JER_EffectiveNP_6__1down",
                          "tree_JET_JER_EffectiveNP_6__1up",
                          "tree_JET_JER_EffectiveNP_7__1down",
                          "tree_JET_JER_EffectiveNP_7__1up",
                          "tree_JET_JER_EffectiveNP_8__1down",
                          "tree_JET_JER_EffectiveNP_8__1up",
                          "tree_JET_JER_EffectiveNP_9__1down",
                          "tree_JET_JER_EffectiveNP_9__1up",
                          "tree_JET_Pileup_OffsetMu__1down",
                          "tree_JET_Pileup_OffsetMu__1up",
                          "tree_JET_Pileup_OffsetNPV__1down",
                          "tree_JET_Pileup_OffsetNPV__1up",
                          "tree_JET_Pileup_PtTerm__1down",
                          "tree_JET_Pileup_PtTerm__1up",
                          "tree_JET_Pileup_RhoTopology__1down",
                          "tree_JET_Pileup_RhoTopology__1up",
                          "tree_JET_PunchThrough_MC16__1down",
                          "tree_JET_PunchThrough_MC16__1up",
                          "tree_JET_SingleParticle_HighPt__1down",
                          "tree_JET_SingleParticle_HighPt__1up",
                          "tree_MET_SoftTrk_ResoPara",
                          "tree_MET_SoftTrk_ResoPerp",
                          "tree_MET_SoftTrk_Scale__1down",
                          "tree_MET_SoftTrk_Scale__1up",
                          "tree_MUON_ID__1down",
                          "tree_MUON_ID__1up",
                          "tree_MUON_MS__1down",
                          "tree_MUON_MS__1up",
                          "tree_MUON_SAGITTA_DATASTAT__1down",
                          "tree_MUON_SAGITTA_DATASTAT__1up",
                          "tree_MUON_SAGITTA_RESBIAS__1down",
                          "tree_MUON_SAGITTA_RESBIAS__1up",
                          "tree_MUON_SCALE__1down",
                          "tree_MUON_SCALE__1up",
                         ]
