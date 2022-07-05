function main_looper(mytree, sumWeight; sfsyst, wgt_factor = 1.0)

    # @hist_prologue
    # @arrow_prologue
    _dict = Dict{Symbol, Hist1D{Float64, Tuple{UnitRange{Int64}}}}()
    for n in (:NN_inZ, :NN_noZ, :NN_DF)
        _dict[Symbol(n, :__NOMINAL)] = Hist1D(Float64; bins=1:2)
        !sfsyst && continue
        for (k,vs) in SF_BRANCH_DICT
            for v in vs, ud in ("1up", "1down")
                _dict[Symbol(n, :__, k, :__, v, :__, ud)] = Hist1D(Float64; bins=1:2)
            end
        end
    end
    hists_dict = dictionary(_dict)

    Threads.@threads for evt in mytree
        ### initial_cut

        v_l_pid = Vcat(evt.v_e_pid, evt.v_m_pid)
        nlepton = length(v_l_pid)
        nlepton != 4 && continue

        
        (; v_e_pt, v_e_eta, v_e_phi, v_e_m,
         v_m_pt, v_m_eta, v_m_phi, v_m_m) = evt
        v_e_tlv = LorentzVectorCyl.(v_e_pt, v_e_eta, v_e_phi, v_e_m)
        v_m_tlv = LorentzVectorCyl.(v_m_pt, v_m_eta, v_m_phi, v_m_m)

        v_l_tlv = Vcat(v_e_tlv, v_m_tlv)

        Z_pair, W_pair, best_Z_mass = Find_Z_Pairs(v_l_pid, v_l_tlv)
        isinf(best_Z_mass) && continue
        wgt = evt.weight / sumWeight * wgt_factor
        other_pair_mass = mass(v_l_tlv[W_pair[1]] + v_l_tlv[W_pair[2]])

        abs(best_Z_mass - Z_m) > 20e3 && continue

        mass_4l = Find_m4l(v_l_tlv)
        mass_4l < 0.0 && continue
        ### end of initial_cut

        !(evt.passTrig) && continue

        v_l_order = sortperm(v_l_tlv; by=pt, rev=true)
        v_l_passIso = get_Isos(evt)

        # `true` in `b_veto` means we've passed the criterial,
        # which means we didn't see a b-tagged
        v_l_medium = Vcat(evt.v_e_LHMedium, evt.v_m_medium) #quality
        v_l_wgtLoose = Vcat(evt.v_e_wgtLoose, evt.v_m_wgtLoose) # quality wgt
        v_l_wgtMedium = Vcat(evt.v_e_wgtMedium, evt.v_m_wgtMedium) #quality wgt

        # W lepton ISO
        v_l_PLTight = Vcat(evt.v_e_passIso_PLImprovedTight, evt.v_m_passIso_PLImprovedTight)
        failed_PLTight = !all(v_l_PLTight[W_pair])
        failed_PLTight && continue
        # v_l_wgtPLTight = Vcat(evt.v_e_wgtIso_PLImprovedTight, evt.v_m_wgtIso_PLImprovedTight)
        # wgt *= reduce(*, v_l_wgtPLTight[W_pair])

        pass_WWZ_cut, wgt, chi2, W_id = WWZ_Cut(
                                                Z_pair,
                                                W_pair,
                                                v_l_pid,
                                                v_l_order,
                                                v_l_wgtLoose,
                                                v_l_wgtMedium,
                                                v_l_tlv,
                                                v_l_passIso,
                                                v_l_medium,
                                                wgt,
                                               )
        !pass_WWZ_cut && continue

        b_wgt, b_veto = Bjet_Cut(evt)
        wgt *= b_wgt
        !b_veto && continue

        l1, l2 = Z_pair
        l3, l4 = W_pair
        v_j_tlv = LorentzVectorCyl.(evt.v_j_pt, evt.v_j_eta, evt.v_j_phi, evt.v_j_m)
        Njet = length(v_j_tlv)

        # Here we calculate the HT's:
        # HT = sum(pt, v_j_tlv; init=0.f0) # hadronic HT
        # leptonic_HT = sum(pt(v_l_tlv[v_l_order[x]]) for x in 1:4)
        # total_HT    = HT + leptonic_HT

        # Distinguish 3 signal channel: inZ, noZ, DF
        if abs(v_l_pid[l3]) != abs(v_l_pid[l4])
            atomic_push!(hists_dict[Symbol(:NN_DF__NOMINAL)], 1, wgt)
        elseif abs(other_pair_mass - Z_m) < 20e3
            atomic_push!(hists_dict[Symbol(:NN_inZ__NOMINAL)], 1, wgt)
        else
            atomic_push!(hists_dict[Symbol(:NN_noZ__NOMINAL)], 1, wgt)
        end

        # @hist_epilogue # -> return hists_dict
        # @arrow_epilogue  # -> return data_ML
    end
    return hists_dict
end
