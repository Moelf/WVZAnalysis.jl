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
