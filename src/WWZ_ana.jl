function WWZ_chi2(v_Z_pair, v_Z_wgt, v_l_pid, v_l_tlv, W_id)
    WWZ_tlv = Vector{LorentzVector{Float64}}(undef, 4)
    pr1 = first(v_Z_pair)
    WWZ_tlv[1] = v_l_tlv[pr1[1]]
    WWZ_tlv[2] = v_l_tlv[pr1[2]]
    WWZ_tlv[3] = v_l_tlv[W_id[1]]
    WWZ_tlv[4] = v_l_tlv[W_id[2]]
    chi2 = 999999.0
    local temp
    Z1_tlv = zero(LorentzVector)
    Z2_tlv = zero(LorentzVector)
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
    v_Z_wgt, v_Z_pair, v_l_pid, v_l_order, v_l_wgt, v_l_tlv, v_l_passIso, v_l_tight, wgt
)
    WWZ_wgt = wgt * first(v_Z_wgt)
    nW = 1
    W_id = Vector{Int}(undef, 2)
    # define W lepton ID and modify weight
    pr1 = first(v_Z_pair)
    @inbounds for vlo in v_l_order
        nW >= 3 && break
        (vlo in pr1) && continue
        WWZ_wgt *= v_l_wgt[vlo]
        W_id[nW] = vlo
        nW += 1
    end
    chi2 = WWZ_chi2(v_Z_pair, v_Z_wgt, v_l_pid, v_l_tlv, W_id)
    FAIL_REUTRN = (false, wgt, Inf, W_id)
    for i in eachindex(v_l_tlv)
        @inbounds for j in (i + 1):length(v_l_tlv)
            v_l_pid[i] + v_l_pid[j] != 0 && continue
            WWZ_dilepton_mass = mass(v_l_tlv[i] + v_l_tlv[j]) / 1000
            WWZ_dilepton_mass < 12 && return FAIL_REUTRN
        end
    end
    pt(v_l_tlv[v_l_order[1]]) < 30e3 && return FAIL_REUTRN
    pt(v_l_tlv[v_l_order[2]]) < 15e3 && return FAIL_REUTRN
    pt(v_l_tlv[v_l_order[3]]) < 8e3 &&  return FAIL_REUTRN
    pt(v_l_tlv[v_l_order[4]]) < 6e3 &&  return FAIL_REUTRN

    # selected lepton min dR
    dR = 999.0
    for i in 1:4
        @inbounds for j in (i + 1):4
            temp = deltaR(v_l_tlv[v_l_order[i]], v_l_tlv[v_l_order[j]])
            dR = ifelse(temp < dR, temp, dR)
        end
    end
    dR < 0.1 && return FAIL_REUTRN

    # summing the charge of the highest pt leptons:
    chargesum = 0
    for i in 1:4
        chargesum += sign(v_l_pid[v_l_order[i]])
    end
    chargesum != 0 && return FAIL_REUTRN

      # Gabriel's best quality (MM)
#     ( (abs(v_l_pid[W_id[1]]) == 11) && !v_l_medium[W_id[1]] ) && return false, wgt, Inf, W_id
#     ( (abs(v_l_pid[W_id[2]]) == 11) && !v_l_medium[W_id[2]] ) && return false, wgt, Inf, W_id
#     ( (abs(v_l_pid[W_id[1]]) == 13) && !v_l_medium[W_id[1]] ) && return false, wgt, Inf, W_id
#     ( (abs(v_l_pid[W_id[2]]) == 13) && !v_l_medium[W_id[2]] ) && return false, wgt, Inf, W_id

      # Gabriel's best isolation (4,1)
#     ( (abs(v_l_pid[W_id[1]]) == 11) && !v_l_passIso[W_id[1]][4] ) && return false, wgt, Inf, W_id
#     ( (abs(v_l_pid[W_id[2]]) == 11) && !v_l_passIso[W_id[2]][4] ) && return false, wgt, Inf, W_id
#     ( (abs(v_l_pid[W_id[1]]) == 13) && !v_l_passIso[W_id[1]][1] ) && return false, wgt, Inf, W_id
#     ( (abs(v_l_pid[W_id[2]]) == 13) && !v_l_passIso[W_id[2]][1] ) && return false, wgt, Inf, W_id

    # following, the quality/isolation cuts from master branch

    #tight cut
    (!v_l_tight[W_id[1]] || !v_l_tight[W_id[2]]) && return FAIL_REUTRN

    #isolation cut, require all 4 of them to be true
    if all((
        v_l_passIso[pr1[1]][1],
        v_l_passIso[pr1[2]][1],
        v_l_passIso[W_id[1]][1],
        v_l_passIso[W_id[2]][1],
    ))
    else
        return FAIL_REUTRN
    end

    return true, WWZ_wgt, chi2, W_id
end # end of WWZ Cut
