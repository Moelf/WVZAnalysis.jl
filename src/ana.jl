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
    return Z_pair, setdiff(1:4, Z_pair), M
end

function Bjet_Cut(evt)
    b = evt.v_j_btag77
    w = evt.v_j_wgt_btag77

    b_wgt = one(eltype(w))
    btag_veto = true
    for (b, w) in zip(b, w)
        b > 0 && (btag_veto = false)
        b_wgt *= w
    end
    return b_wgt, btag_veto
end

function main_looper(r::ROOTFile, sumWeight; shape_variation = "NOMINAL", kw...)
    mytree = LazyTree(r, "tree_" * shape_variation)
    return main_looper(mytree, sumWeight; shape_variation, kw...)
end
