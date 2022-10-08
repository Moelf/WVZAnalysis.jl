function ZZZ_Cut(v_Z_pair, v_ignore, v_l_pid, v_l_tlv, wgt)
    nlepton = length(v_l_pid)
    nZ = length(v_Z_pair)
    nlepton < 6 && return false, wgt

    nZ == 1 && return false, wgt

    if nZ == 2
        # recovery? doesn't do anything other than update cut flow
    end

    nlepton > 6 && return false, wgt

    # 3 SFOS
    if nZ == 3
        # update weights?
    end

    if nZ == 2
        Z3_tlv = zero(LorentzVector)
        for i in eachindex(v_l_pid)
            i ∈ v_ignore && continue
            Z3_tlv += v_l_tlv[i]
        end
        m_Z3 = mass(Z3_tlv)

        m_Z3 <= 40e3 && return false, wgt
    end

    if nZ == 3
        sp = v_Z_pair[3]
        m_2ndpair = mass(v_l_tlv[sp[1]] + v_l_tlv[sp[2]])
        m_2ndpair <= 40e3 && return false, wgt
    end

    return true, wgt
end # end of ZZZ_Cut