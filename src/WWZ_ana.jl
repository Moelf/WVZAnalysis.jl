function WWZ_chi2(pr1, W_id, v_l_pid, v_l_tlv)
    WWZ_tlv = Vector{eltype(v_l_tlv)}(undef, 4)
    WWZ_tlv[1] = v_l_tlv[pr1[1]]
    WWZ_tlv[2] = v_l_tlv[pr1[2]]
    WWZ_tlv[3] = v_l_tlv[W_id[1]]
    WWZ_tlv[4] = v_l_tlv[W_id[2]]
    chi2 = 999999.0
    local temp
    Z1_tlv = zero(WWZ_tlv[1])
    Z2_tlv = zero(WWZ_tlv[1])
    for i in 2:4
        Z1_tlv = WWZ_tlv[1] + WWZ_tlv[i]
        for j in 2:4
            (j == i) && continue
            Z2_tlv += WWZ_tlv[j]
        end
        mz1 = mass(Z1_tlv)
        mz2 = mass(Z2_tlv)
        temp =
            (
                (mz1 - Z_m)^2 + (mz2 - Z_m)^2
            ) / 2495.2^2
        (temp < chi2) && (chi2 = temp)
    end
    return chi2
end

Base.@propagate_inbounds function WWZ_Cut(
    Z_pair, W_pair, v_l_pid, v_l_order, v_l_wgtLoose, v_l_wgtMedium, v_l_tlv, v_l_passIso, v_l_medium, wgt
)
    nW = 1
    chi2 = Inf
    # chi2 = WWZ_chi2(Z_pair, W_pair, v_l_pid, v_l_tlv)
    FAIL = (false, wgt, Inf, W_pair)
    for i in eachindex(v_l_tlv)
        @inbounds for j in (i + 1):length(v_l_tlv)
            v_l_pid[i] + v_l_pid[j] != 0 && continue
            WWZ_dilepton_mass = mass(v_l_tlv[i] + v_l_tlv[j]) / 1000
            WWZ_dilepton_mass < 12 && return FAIL
        end
    end
    pt(v_l_tlv[v_l_order[1]]) < 30e3 && return FAIL
    pt(v_l_tlv[v_l_order[2]]) < 15e3 && return FAIL
    pt(v_l_tlv[v_l_order[3]]) < 8e3 &&  return FAIL
    pt(v_l_tlv[v_l_order[4]]) < 6e3 &&  return FAIL

    # selected lepton min dR
    dR = 999.0
    for i in 1:4
        @inbounds for j in (i + 1):4
            temp = deltaR(v_l_tlv[v_l_order[i]], v_l_tlv[v_l_order[j]])
            dR = ifelse(temp < dR, temp, dR)
        end
    end
    dR < 0.1 && return FAIL

    # summing the charge of the highest pt leptons:
    chargesum = 0
    for i in 1:4
        chargesum += sign(v_l_pid[v_l_order[i]])
    end
    chargesum != 0 && return FAIL

    for i in 1:2
        ### for Z bosons
        # Overall best quality (Loose(e) and Loose(mu))
        # Isolation (5,3)
        ( (abs(v_l_pid[Z_Pair[i]]) == 11) && !v_l_passIso[Z_Pair[i]][3] ) && return FAIL
        ( (abs(v_l_pid[Z_Pair[i]]) == 13) && !v_l_passIso[Z_Pair[i]][2] ) && return FAIL
        ### for W bosons
        # Overall best quality (Medium(e) and Medium(mu))
        ( !v_l_medium[W_pair[i]] ) && return FAIL
    end
    # define W lepton ID and modify weight
    WWZ_wgt = wgt * v_l_wgtLoose[Z_pair[1]] * v_l_wgtLoose[Z_pair[2]] * v_l_wgtLoose[W_pair[1]] * v_l_wgtLoose[W_pair[2]]
    return true, WWZ_wgt, chi2, W_pair
end # end of WWZ Cut
