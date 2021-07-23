const Z_m = 91.1876 * 10^3 # in MeV

mass(lv::LorentzVector) = sqrt(dot(lv, lv))
pt(lv::LorentzVector) = sqrt(lv.x^2 + lv.y^2)
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
    fZ == 0.0 && return 0.0
    # Warning("PseudoRapidity","transvers momentum = 0! return +/- 10e10");
    fZ > 0.0 && return 10e10
    return -10e10
end
function phi(lv::LorentzVector)
    return (lv.x == 0.0 && lv.y == 0.0) ? 0.0 : atan(lv.y, lv.x)
end

function phi_mpi_pi(x)
    while (x >= pi)
        x -= 2pi
    end
    while (x < -pi)
        x += 2pi
    end
    return x
end

function deltaR(lv1::LorentzVector, lv2::LorentzVector)
    deta = eta(lv1) - eta(lv2)
    dphi = phi_mpi_pi(phi(lv1) - phi(lv2))
    return sqrt(deta * deta + dphi * dphi)
end
