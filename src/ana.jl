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
    main_looper(mytree, sumWeight)
end

function main_looper(mytree, sumWeight)
    hists_dict = Dict{Symbol,Hist1D}(
        :Z_mass_first => Hist1D(Float32; bins=0:10:200),
        :WZZ_ZZ_mass => Hist1D(Float32; bins=0:10:800),
        :WWZ_MET => Hist1D(Float32; bins=0:5:400),
    )
    @inbounds for (i, evt) in enumerate(mytree)
        ### initial_cut
        e_mask = evt.v_e_fwd
        e_mask .⊻= true
        m_mask = evt.v_m_lowpt
        m_mask .⊻= true

        v_l_pid = vcat(evt.v_e_pid[e_mask], evt.v_m_pid[m_mask])
        nlepton = length(v_l_pid)
        nlepton <= 3 && continue

        v_l_tlv = vcat(evt.v_e_tlv[e_mask], evt.v_m_tlv[m_mask])
        v_l_wgt = vcat(evt.v_e_wgtLoose[e_mask], evt.v_m_wgtLoose[m_mask])

        v_Z_pair, v_Z_wgt, v_ignore = Find_Z_Pairs(v_l_pid, v_l_tlv, v_l_wgt)
        isempty(v_Z_pair) && continue
        zpr1 = first(v_Z_pair)
        wgt = evt.weight / sumWeight
        push!(
            hists_dict[:Z_mass_first], mass(v_l_tlv[zpr1[1]] + v_l_tlv[zpr1[2]]) / 1000, wgt
        )

        abs(mass(v_l_tlv[zpr1[1]] + v_l_tlv[zpr1[2]]) - Z_m) > 20e3 && continue

        v_l_order = sortperm(v_l_tlv; by=pt, rev=true)
        mass_4l = Find_m4l(v_Z_pair, v_l_tlv, v_l_order)
        mass_4l < 0.0 && continue
        ### end of initial_cut

        pass_ZZZ_cut, wgt = ZZZ_Cut(v_Z_pair, v_ignore, v_l_pid, v_l_tlv, wgt)
        if pass_ZZZ_cut
            continue
        end
        !(evt.passTrig) && continue

        v_l_passIso = get_Isos(e_mask, m_mask, evt)

        pass_WZZ_cut, wgt = WZZ_Cut(
            v_Z_wgt, v_Z_pair, v_l_pid, v_l_order, v_l_wgt, v_l_tlv, v_l_passIso, wgt
        )
        # `true` in `b_veto` means we've passed the criterial,
        # which means we didn't see a b-tagged
        if pass_WZZ_cut
            zpr2 = v_Z_pair[2]
            Z1_vlt = v_l_tlv[zpr1[1]] + v_l_tlv[zpr1[2]]
            Z2_vlt = v_l_tlv[zpr2[1]] + v_l_tlv[zpr2[2]]
            ZZ_vlt = Z1_vlt + Z2_vlt
            ZZ_mass = mass(ZZ_vlt)

            push!(hists_dict[:WZZ_ZZ_mass], ZZ_mass / 1000, wgt)
            continue
        end

        v_l_tight = vcat(evt.v_e_LHTight[e_mask], evt.v_m_tight[m_mask])
        pass_WWZ_cut, wgt = WWZ_Cut(
            v_Z_wgt,
            v_Z_pair,
            v_l_pid,
            v_l_order,
            v_l_wgt,
            v_l_tlv,
            v_l_passIso,
            v_l_tight,
            wgt,
        )
        b_wgt, b_veto = Bjet_Cut(evt)
        if pass_WWZ_cut
            if (b_veto)
                wgt *= b_wgt
            else
                continue
            end
            # @show i
            # @show wgt
            push!(hists_dict[:WWZ_MET], evt.MET / 1000, wgt)
            continue
        end
    end
    return hists_dict
end
