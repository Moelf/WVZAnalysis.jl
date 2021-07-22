function WZZ_Cut(nlepton, nZ, v_Z_wgt, v_Z_pair, v_l_pid, v_l_order, v_l_wgt, v_l_tlv, wgt)
    nlepton < 5 && return false

    nZ < 2 && return false

    WZZ_wgt = wgt * v_Z_wgt[1] * v_Z_wgt[2]

    pr1 = v_Z_pair[1]
    pr2 = v_Z_pair[2]
    for i in eachindex(v_l_pid)
        vlo = v_l_order[i]
        if (vlo == pr1[1] || vlo == pr1[2] || vlo == pr2[1] || vlo == pr2[2])
            continue
        end
        WZZ_wgt *= v_l_wgt[vlo]
        break
    end

    abs(mass((v_l_tlv[pr2[1]]+v_l_tlv[pr2[2]])) - Z_m) > 20e3 && return false;

    pt(v_l_tlv[v_l_order[1]]) < 50e3 && return false
    pt(v_l_tlv[v_l_order[2]]) < 40e3 && return false
    pt(v_l_tlv[v_l_order[3]]) < 20e3 && return false

    # FIXME
    # Bjet_Cut("WZZ","WZZ",WZZ_wgt);

    return true
end # end of WZZ Cut
