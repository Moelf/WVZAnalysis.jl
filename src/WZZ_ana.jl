function WZZ_Cut(v_Z_wgt, v_Z_pair, v_l_pid, v_l_order, v_l_wgt, v_l_tlv, wgt)
    nlepton = length(v_l_pid)
    fifth_l = -1
    nZ = length(v_Z_pair)

    nlepton < 5 && return false, wgt

    nZ < 2 && return false, wgt

    WZZ_wgt = wgt * v_Z_wgt[1] * v_Z_wgt[2]

    pr1 = v_Z_pair[1]
    pr2 = v_Z_pair[2]
    for vlo in v_l_order
        if (vlo == pr1[1] || vlo == pr1[2] || vlo == pr2[1] || vlo == pr2[2])
            continue
        end
        fifth_l = vlo
        WZZ_wgt *= v_l_wgt[vlo]
        break
    end

    abs(mass((v_l_tlv[pr2[1]]+v_l_tlv[pr2[2]])) - Z_m) > 20e3 && return false, wgt

    return true, WZZ_wgt
end # end of WZZ Cut
