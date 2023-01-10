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
