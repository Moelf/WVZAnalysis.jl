const Z_m = 91.1876 * 10^3 # in MeV

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

function deltaR(lv1::LorentzVector, lv2::LorentzVector)
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
