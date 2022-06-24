function main_looper(mytree, sumWeight; sfsyst, wgt_factor = 1.0)

    #println("Running mainlooper")
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
    MET_hist = Hist1D(Float64; bins=0:4:150)
    Z_mass_hist = Hist1D(Float64; bins=50:2:150)
    bad_mass_hist = Hist1D(Float64; bins=50:4:150)
    eta_hist = Hist1D(Float64; bins=-4:0.2:4)
    phi_hist = Hist1D(Float64; bins=-4:0.2:4)
    pt_hist = Hist1D(Float64; bins=0:4:200)


    Threads.@threads for evt in mytree
        #println("v_e_pid: ",(evt.v_e_pid))
        ### initial_cut
        v_l_pid = Vcat(evt.v_e_pid, evt.v_m_pid)
        nlepton = length(v_l_pid)
        nlepton != 4 && continue

        # (;v_e_tlv, v_m_tlv) = evt
        
        (; v_e_pt, v_e_eta, v_e_phi, v_e_m,
         v_m_pt, v_m_eta, v_m_phi, v_m_m) = evt
        v_e_tlv = LorentzVectorCyl.(v_e_pt, v_e_eta, v_e_phi, v_e_m)
        v_m_tlv = LorentzVectorCyl.(v_m_pt, v_m_eta, v_m_phi, v_m_m)

        v_l_tlv = Vcat(v_e_tlv, v_m_tlv)
        v_l_eta = Vcat(v_e_eta,v_m_eta)
        v_l_phi = Vcat(v_e_phi,v_m_phi)
        v_l_pt = Vcat(v_e_pt,v_m_pt)
        v_l_wgt = Vcat(evt.v_e_wgtLoose, evt.v_m_wgtLoose)

        zpr1, other, best_Z_mass = Find_Z_Pairs(v_l_pid, v_l_tlv)
        isinf(best_Z_mass) && continue
        wgt = evt.weight / sumWeight * wgt_factor

	other_pair_mass = mass(v_l_tlv[other[1]] + v_l_tlv[other[2]]) 
    
        mass_4l = Find_m4l(v_l_tlv)
        mass_4l < 0.0 && continue
        ### end of initial_cut

        20e3 > abs(best_Z_mass - Z_m) && continue
        40e3 < abs(best_Z_mass - Z_m) && continue
        !(evt.passTrig) && continue
       
        v_l_order = sortperm(v_l_tlv; by=pt, rev=true)
        v_l_passIso = get_Isos(evt)
        
        # `true` in `b_veto` means we've passed the criterial,
        # which means we didn't see a b-tagged
        v_l_medium = Vcat(evt.v_e_LHMedium, evt.v_m_medium)
        pass_WWZ_cut, wgt, chi2, W_id = WWZ_Cut(
                                                zpr1,
                                                other,
                                                v_l_pid,
                                                v_l_order,
                                                v_l_wgt,
                                                v_l_tlv,
                                                v_l_passIso,
                                                v_l_medium,
                                                wgt,
                                               )
        b_wgt, b_veto = Bjet_Cut(evt)
        !pass_WWZ_cut && continue
        !b_veto && continue
        wgt *= b_wgt
        l1, l2 = zpr1
        l3, l4 = W_id
        # v_j_tlv = evt.v_j_tlv
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

        atomic_push!(MET_hist, evt.MET/1000, wgt)
        atomic_push!(Z_mass_hist,best_Z_mass/1000, wgt)
	atomic_push!(eta_hist,eta(v_l_tlv[zpr1[1]]+v_l_tlv[zpr1[2]]), wgt)
        atomic_push!(phi_hist,phi(v_l_tlv[zpr1[1]]+v_l_tlv[zpr1[2]]), wgt)
        atomic_push!(pt_hist,pt(v_l_tlv[zpr1[1]]+v_l_tlv[zpr1[2]])/1000, wgt)
        atomic_push!(bad_mass_hist,other_pair_mass/1000, wgt)
    end
    return MET_hist, Z_mass_hist, eta_hist, phi_hist, pt_hist, bad_mass_hist
    #return len(eta)
end

