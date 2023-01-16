function Find_Z_Pairs(v_l_pids, v_l_tlv)
    Z_pair = [-1,-1]
    idxs = eachindex(v_l_pids)
    M = Inf
    for i in idxs
        vi = v_l_tlv[i]
        pidi = v_l_pids[i]
        for j in (i + 1):lastindex(v_l_pids)
            pidi != -v_l_pids[j] && continue # require OS
            m0 = mass(vi + v_l_tlv[j])
            if abs(m0 - Z_m) < abs(M - Z_m)
                M = m0
                Z_pair[1] = i
                Z_pair[2] = j
            end
        end
    end
    W_pair = setdiff(1:4, Z_pair)

    # order W lep by pT
    if pt(v_l_tlv[W_pair[1]]) < pt(v_l_tlv[W_pair[2]])
        W_pair[1], W_pair[2] = W_pair[2], W_pair[1]
    end

    return Z_pair, W_pair, M
end

function WWZ_chi2(Z_pair, W_pair, v_l_pid, v_l_tlv)
    WWZ_tlv = Vector{eltype(v_l_tlv)}(undef, 4)
    WWZ_tlv[1] = v_l_tlv[Z_pair[1]]
    WWZ_tlv[2] = v_l_tlv[Z_pair[2]]
    WWZ_tlv[3] = v_l_tlv[W_pair[1]]
    WWZ_tlv[4] = v_l_tlv[W_pair[2]]
    chi2 = Inf
    local temp
    Z1_tlv = zero(WWZ_tlv[1])
    Z2_tlv = zero(WWZ_tlv[1])
    @inbounds for i in 2:4
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
