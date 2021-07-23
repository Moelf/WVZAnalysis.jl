function Find_Z_Pairs(v_l_pids, v_l_tlv, v_l_wgt)
    v_Z_pair = Tuple{Int,Int}[]
    v_Z_wgt = Float32[]

    v_ignore = Set{Int}()
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
        push!(v_ignore, temp_tup[1])
        push!(v_ignore, temp_tup[2])

        push!(v_Z_pair, temp_tup)
        push!(v_Z_wgt, v_l_wgt[temp_tup[1]] * v_l_wgt[temp_tup[2]])
    end

    return v_Z_pair, v_Z_wgt, v_ignore
end

function Find_m4l(v_Z_pair, v_l_tlv, v_l_order)
    m4l = first(v_Z_pair)

    @inbounds for i in eachindex(v_l_order)
        length(m4l) >= 4 && break
        (v_l_order[i] == m4l[1] || v_l_order[i] == m4l[2]) && continue
        m4l = (m4l..., v_l_order[i])
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
    tagnums = ("60", "70", "77", "85")
    _btags = [getproperty(evt, Symbol(:v_j_btag, i)) for i in tagnums]
    _wgts = [getproperty(evt, Symbol(:v_j_wgt_btag, i)) for i in tagnums]

    b_wgt = Vector{Float64}(undef, 4)
    b_veto = Vector{Bool}(undef, 4)

    for idx in eachindex(_btags, _wgts)
        btag_wgt = 1.0
        btag_veto = true
        for (b, w) in zip(_btags[idx], _wgts[idx])
            b > 0 && (btag_veto = false)
            btag_wgt *= w
        end
        # normally there's a cutflow line here
        b_wgt[idx] = btag_wgt
        b_veto[idx] = btag_veto
    end
    return b_wgt, b_veto
end

function main_looper(r::ROOTFile)
    sumWeight = r["sumWeight"][:fN][3]
    mytree = LazyTree(
        r,
        "tree_NOMINAL",
        ["passTrig", r"v_j_(wgt_)?btag.*", "weight", r"v_(e|m|j)_(fwd|tlv|wgtLoose|pid|lowpt)$"],
    )
    hist_WZZ_ZZ_mass = Hist1D(; bins=0:10:800)
    hist_Z_mass_first = Hist1D(Float32; bins=0:10:200)
    for evt in mytree
        ### initial_cut
        wgt = evt.weight / sumWeight
        e_mask = .!(evt.v_e_fwd)
        m_mask = .!(evt.v_m_lowpt)
        v_l_pid = vcat(evt.v_e_pid[e_mask], evt.v_m_pid[m_mask])
        nlepton = length(v_l_pid)
        nlepton <= 3 && continue

        v_l_tlv = vcat(evt.v_e_tlv[e_mask], evt.v_m_tlv[m_mask])
        v_l_wgt = vcat(evt.v_e_wgtLoose[e_mask], evt.v_m_wgtLoose[m_mask])

        v_Z_pair, v_Z_wgt, v_ignore = Find_Z_Pairs(v_l_pid, v_l_tlv, v_l_wgt)
        isempty(v_Z_pair) && continue
        zpr1 = first(v_Z_pair)
        push!(hist_Z_mass_first, mass(v_l_tlv[zpr1[1]] + v_l_tlv[zpr1[2]])/1000 ,wgt)

        best_p = first(v_Z_pair)
        abs(mass(v_l_tlv[best_p[1]] + v_l_tlv[best_p[2]]) - Z_m) > 20e3 && continue

        v_l_order = sortperm(v_l_tlv; by=mass)
        mass_4l = Find_m4l(v_Z_pair, v_l_tlv, v_l_order)
        # end of initial_cut

        pass_ZZZ_cut, wgt = ZZZ_Cut(v_Z_pair, v_ignore, v_l_pid, v_l_tlv, wgt)
        if pass_ZZZ_cut
            continue
        end
        !(evt.passTrig) && continue

        pass_WZZ_cut, wgt = WZZ_Cut(
            v_Z_wgt, v_Z_pair, v_l_pid, v_l_order, v_l_wgt, v_l_tlv, wgt
        )
        # `true` in `b_veto` means we've passed the criterial,
        # which means we didn't see a b-tagged
        b_wgt, b_veto = Bjet_Cut(evt)
        if (b_veto[3])
            wgt *= b_wgt[3]
        else
            continue
        end
        if pass_WZZ_cut
            zpr2 = v_Z_pair[2]
            Z1_vlt = v_l_tlv[zpr1[1]] + v_l_tlv[zpr1[2]]
            Z2_vlt = v_l_tlv[zpr2[1]] + v_l_tlv[zpr2[2]]
            ZZ_vlt = Z1_vlt + Z2_vlt
            ZZ_mass = mass(ZZ_vlt)
            ZZ_pt = pt(ZZ_vlt)

            # @show ZZ_mass / 1000
            push!(hist_WZZ_ZZ_mass, ZZ_mass / 1000, wgt)
            continue
        end

        pass_WWZ_cut, wgt = WWZ_Cut(
            v_Z_wgt, v_Z_pair, v_l_pid, v_l_order, v_l_wgt, v_l_tlv, wgt
        )
        if pass_WWZ_cut
            continue
        end

    end
    return hist_WZZ_ZZ_mass
end
