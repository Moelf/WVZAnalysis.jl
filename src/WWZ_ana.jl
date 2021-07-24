function WWZ_chi2(v_Z_pair, v_Z_wgt, v_l_pid, v_l_tlv, W_id)
    WWZ_tlv = Vector{LorenztVector{Float64}}(undef, 4)
    pr1 = first(v_Z_pair)
    WWZ_tlv[1] = v_l_tlv[pr1[1]]
    WWZ_tlv[2] = v_l_tlv[pr1[2]]
    WWZ_tlv[3] = v_l_tlv[W_id[1]]
    WWZ_tlv[4] = v_l_tlv[W_id[2]]
    chi2 = 999999.0
    local temp
    local Z1_tlv , Z2_tlv
    for i in 2:4
        Z1_tlv = WWZ_tlv[1] + WWZ_tlv[i]
        for j in 2:4
            (j == i) && continue
            Z2_tlv += WWZ_tlv[j]
        end
        temp =
            (
                (mass(Z1_tlv) - Z_m) * (mass(Z1_tlv) - Z_m) +
                (mass(Z2_tlv) - Z_m) * (mass(Z2_tlv) - Z_m)
            ) / (2495.2 * 2495.2)
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
    # chi2 = WWZ_chi2(v_Z_pair, v_Z_wgt, v_l_pid, v_l_tlv, W_id)
    for i in eachindex(v_l_tlv)
        @inbounds for j in (i + 1):length(v_l_tlv)
            v_l_pid[i] + v_l_pid[j] != 0 && continue
            WWZ_dilepton_mass = mass(v_l_tlv[i] + v_l_tlv[j]) / 1000
            WWZ_dilepton_mass < 12 && return false, wgt
        end
    end
    pt(v_l_tlv[v_l_order[1]]) < 30e3 && return false, wgt
    pt(v_l_tlv[v_l_order[2]]) < 15e3 && return false, wgt
    pt(v_l_tlv[v_l_order[3]]) < 8e3 && return false, wgt
    pt(v_l_tlv[v_l_order[4]]) < 6e3 && return false, wgt

    # selected lepton min dR
    dR = 999.0
    for i in 1:4
        @inbounds for j in (i + 1):4
            temp = deltaR(v_l_tlv[v_l_order[i]], v_l_tlv[v_l_order[j]])
            dR = ifelse(temp < dR, temp, dR)
        end
    end
    dR < 0.1 && return false, wgt

    #tight cut
    (!v_l_tight[W_id[1]] || !v_l_tight[W_id[2]]) && return false, wgt

    #isolation cut, require all 4 of them to be true
    if all((
        v_l_passIso[pr1[1]][1],
        v_l_passIso[pr1[2]][1],
        v_l_passIso[W_id[1]][1],
        v_l_passIso[W_id[2]][1],
    ))
    else
        return false, wgt
    end

    return true, WWZ_wgt
end # end of WWZ Cut
