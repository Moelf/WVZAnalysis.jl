function Find_Z_Pairs(v_l_pids, v_l_tlv, v_l_wgt)
    v_Z_pair = Tuple{Int,Int}[]
    v_Z_wgt = Float32[]

    v_ignore = Int[]
    @inbounds while length(v_ignore) < length(v_l_pids)
        M = Inf
        local temp_tup
        for i in eachindex(v_l_pids) # electrons loop
            i ∈ v_ignore && continue
            for j in (i + 1):length(v_l_pids)
                j ∈ v_ignore && continue
                v_l_pids[i] != -v_l_pids[j] && continue # require OS
                m0 = mass(v_l_tlv[i] + v_l_tlv[j])
                if abs(m0 - Z_m) < abs(M - Z_m)
                    M = m0
                    temp_tup = (i, j)
                end
            end
        end
        isinf(M) && break # can't find any more pairs
        push!(v_ignore, temp_tup...)

        push!(v_Z_pair, temp_tup)
        push!(v_Z_wgt, v_l_wgt[temp_tup[1]] * v_l_wgt[temp_tup[2]])
    end

    return v_Z_pair, v_Z_wgt, v_ignore
end

function Find_m4l(v_Z_pair, v_l_tlv, v_l_order)
    m4l = first(v_Z_pair)

    @inbounds for vlo in v_l_order
        length(m4l) >= 4 && break
        (vlo ∈ m4l) && continue
        m4l = (m4l..., vlo)
    end
    # require 2SFOS
    length(m4l) != 4 && return -1.0

    tlv_4l = zero(LorentzVector)
    for idx in m4l
        tlv_4l += v_l_tlv[idx]
    end

    return mass(tlv_4l)
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

main_looper(s::String) = main_looper(ROOTFile(s))

function main_looper(r::ROOTFile)
    sumWeight = r["sumWeight"][:fN][3]
    mytree = LazyTree(
        r,
        "tree_NOMINAL",
        [
            "MET",
            "passTrig",
            r"v_(e|m)_(LHTight|tight)",
            r"v_j_(wgt_)?btag.*",
            r"v_(e|m)_passIso_.*",
            "weight",
            r"v_(e|m|j)_(fwd|tlv|wgtLoose|pid|lowpt)$",
        ],
    )
    return main_looper(mytree, sumWeight)
end

