function ZZZ_Cut(v_Z_pair, nZ, nlepton, v_ignore, v_l_pid, v_l_tlv)
    nlepton < 6 && return false

    nZ == 1 && return false

    if nZ == 2
        # recovery? doesn't do anything other than update cut flow
    end

    nlepton > 6 && return false

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

        m_Z3 < 40e3 && return false
    end

    if nZ == 3
        sp = v_Z_pair[2]
        m_2ndpair = mass(v_l_tlv[sp[1]] + v_l_tlv[sp[2]])
        m_2ndpair <= 40e3 && return false
    end

    return true
end # end of ZZZ_Cut
