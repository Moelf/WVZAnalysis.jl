function main_looper(mytree, sumWeight; sfsyst, wgt_factor = 1.0, arrow_making=false, isdata=false, controlregion=nothing)
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
    dict, executor = if arrow_making
        arrow_init(), NonThreadedEx()
    else
        kinematic_hist_init(), ThreadedEx()
    end

    model, scales, minimums = init_ONNX()

    @floop executor for evt in mytree
    # for evt in mytree
        ### initial_cut
        v_m_eta_orig, v_e_caloeta_orig = evt.v_m_eta, evt.v_e_caloeta
        e_etamask = [abs(η) < 2.47 && (abs(η)<1.37 || abs(η)>1.52) for η in v_e_caloeta_orig]
        m_etamask = [abs(η) < 2.5 for η in v_m_eta_orig]
        v_l_pid = @views Vcat(evt.v_e_pid[e_etamask], evt.v_m_pid[m_etamask])
        Nlep = length(v_l_pid)
        Nlep != 4 && continue
        (; v_e_pt, v_e_eta, v_e_phi, v_e_m,
         v_m_pt, v_m_eta, v_m_phi, v_m_m) = evt
        v_e_tlv = @views LorentzVectorCyl.(v_e_pt[e_etamask], v_e_eta[e_etamask], v_e_phi[e_etamask], v_e_m[e_etamask])
        v_m_tlv = @views LorentzVectorCyl.(v_m_pt[m_etamask], v_m_eta[m_etamask], v_m_phi[m_etamask], v_m_m[m_etamask])
        v_l_tlv = Vcat(v_e_tlv, v_m_tlv)
        v_l_eta = Vcat(v_e_eta,v_m_eta)
        v_l_phi = Vcat(v_e_phi,v_m_phi)
        v_l_pt = Vcat(v_e_pt,v_m_pt)
        v_l_wgt = Vcat(evt.v_e_wgtLoose, evt.v_m_wgtLoose)
        Z_pair, W_pair, best_Z_mass = Find_Z_Pairs(v_l_pid, v_l_tlv)
        isinf(best_Z_mass) && continue
        other_mass = mass(v_l_tlv[W_pair[1]] + v_l_tlv[W_pair[2]])
        wgt = evt.weight / sumWeight * wgt_factor
        abs(best_Z_mass - Z_m) > 20e3 && continue
        mass_4l = mass(sum(v_l_tlv))
        mass_4l < 0.0 && continue
        ### end of initial_cut
        !(evt.passTrig) && continue
        v_l_order = sortperm(v_l_tlv; by=pt, rev=true)
        wgt = evt.weight / sumWeight * wgt_factor
        v_l_medium =    @views Vcat(evt.v_e_LHMedium[e_etamask] , evt.v_m_medium[m_etamask]) #quality
        v_l_wgtLoose =  @views Vcat(evt.v_e_wgtLoose[e_etamask] , evt.v_m_wgtLoose[m_etamask]) # quality wgt
        v_l_wgtMedium = @views Vcat(evt.v_e_wgtMedium[e_etamask], evt.v_m_wgtMedium[m_etamask]) #quality wgt
        ############## use PLIV for W lepton ISO #################
        v_l_PLTight = Vcat(evt.v_e_passIso_PLImprovedTight[e_etamask], evt.v_m_passIso_PLImprovedTight[m_etamask])
        if controlregion == :Zjets
            failed_PLTight = any(==(true), v_l_PLTight[W_pair])
            failed_PLTight && continue
            v_l_Loose = Vcat(evt.v_e_passIso_Loose_VarRad, evt.v_m_passIso_PflowLoose_VarRad)
            failed_Loose = all(==(true), v_l_Loose[W_pair])
            failed_Loose && continue
        else
            failed_PLTight = any(==(false), v_l_PLTight[W_pair])
            failed_PLTight && continue
            v_l_wgtPLTight = Vcat(evt.v_e_wgtIso_PLImprovedTight_Medium[e_etamask], evt.v_m_wgtIso_PLImprovedTight[m_etamask])
            wgt *= @views reduce(*, v_l_wgtPLTight[W_pair]; init = 1.0)
        end
        wgt *= @views reduce(*, evt.v_e_wgtReco[e_etamask])
        wgt *= @views reduce(*, evt.v_m_wgtTTVA[m_etamask])
        (;
         v_e_passIso_Loose_VarRad,
         v_m_passIso_PflowLoose_VarRad,
         v_e_wgtIso_Loose_VarRad_LooseBLayer,
         v_m_wgtIso_PflowLoose_VarRad,
        ) = evt
        v_l_passIso = @views Vcat(v_e_passIso_Loose_VarRad[e_etamask], v_m_passIso_PflowLoose_VarRad[m_etamask])
        v_l_wgtIso =  @views Vcat(v_e_wgtIso_Loose_VarRad_LooseBLayer[e_etamask], v_m_wgtIso_PflowLoose_VarRad[m_etamask])
        pass_WWZ_cut, wgt, chisq, W_id = WWZ_Cut(
                                                 Z_pair,
                                                 W_pair,
                                                 v_l_pid,
                                                 v_l_order,
                                                 v_l_wgtLoose,
                                                 v_l_medium,
                                                 v_l_wgtMedium,
                                                 v_l_tlv,
                                                 v_l_passIso,
                                                 v_l_wgtIso,
                                                 wgt,
                                                 isdata
                                                )
        !pass_WWZ_cut && continue
        # in `b_veto == true` means we've passed the criterial,
        # which means we didn't see a b-jet
        b_wgt, b_veto = Bjet_Cut(evt)
        (; MET, METSig, METPhi) = evt
        MET /= 1000
        if controlregion == :ttZ
            MET < 20 && continue
            (other_mass > 80e3 && other_mass < 100e3) && continue
            b_veto && continue
        else
            !b_veto && continue
        end
        wgt *= b_wgt
        # Distinguish 3 signal channel: inZ, noZ, DF
        l3, l4 = W_pair
        SR = if abs(v_l_pid[l3]) != abs(v_l_pid[l4]) #DF
            2
        elseif abs(other_mass - Z_m) < 20e3 # SF_inZ
            0
        else #SF_noZ
            1
        end
        # force wgt to 1 for data
        if isdata
            wgt = 1.0
        end
        lep1_pid, lep2_pid, lep3_pid, lep4_pid = @view v_l_pid[v_l_order]
        pt_1, pt_2, pt_3, pt_4 = pt.(@view v_l_tlv[v_l_order]) ./ 1000
        eta_1, eta_2, eta_3, eta_4 = eta.(@view v_l_tlv[v_l_order])
        phi_1, phi_2, phi_3, phi_4 = phi.(@view v_l_tlv[v_l_order])
        phi_1, phi_2, phi_3, phi_4 = phi.(@view v_l_tlv[v_l_order])
        v_j_tlv = LorentzVectorCyl.(evt.v_j_pt, evt.v_j_eta, evt.v_j_phi, evt.v_j_m)
        v_j_order = sortperm(v_j_tlv; by=pt, rev=true)
        Njet = length(v_j_tlv)
        Zcand_mass = best_Z_mass / 1000
        Z_eta = eta(v_l_tlv[Z_pair[1]]+v_l_tlv[Z_pair[2]])
        Z_phi = phi(v_l_tlv[Z_pair[1]]+v_l_tlv[Z_pair[2]])
        Z_pt = pt(v_l_tlv[Z_pair[1]]+v_l_tlv[Z_pair[2]])/1000
        HT = sum(pt, v_j_tlv; init=0.f0)/1000 # hadronic HT
        leptonic_HT = sum(pt(v_l_tlv[v_l_order[x]]) for x in 1:4)/1000
        total_HT    = HT + leptonic_HT
        total_events = 1
  
        Z_tlv = v_l_tlv[Z_pair[1]]+v_l_tlv[Z_pair[2]]
        Z_rapidity = 0.5 * log((energy(Z_tlv)+pz(Z_tlv))/(energy(Z_tlv)-pz(Z_tlv)))
        other_mass /= 1000
        mass_4l /= 1000
        Zlep1_pt, Zlep2_pt = pt.(@view v_l_tlv[Z_pair]) ./ 1000
        Zlep1_eta, Zlep2_eta = eta.(@view v_l_tlv[Z_pair])
        Zlep1_phi, Zlep2_phi = phi.(@view v_l_tlv[Z_pair])
        Zlep1_pid, Zlep2_pid = @view v_l_pid[Z_pair]
        Zlep1_dphi = LorentzVectorHEP.phi_mpi_pi(Z_phi - phi(v_l_tlv[Z_pair[1]]))
        Zlep2_dphi = LorentzVectorHEP.phi_mpi_pi(Z_phi - phi(v_l_tlv[Z_pair[2]]))
        Wlep1_pt, Wlep2_pt = pt.(@view v_l_tlv[W_pair]) ./ 1000
        Wlep1_eta, Wlep2_eta = eta.(@view v_l_tlv[W_pair])
        Wlep1_phi, Wlep2_phi = phi.(@view v_l_tlv[W_pair])
        Wlep1_pid, Wlep2_pid = @view v_l_pid[W_pair]
        Wlep1_dphi = LorentzVectorHEP.phi_mpi_pi(Z_phi - phi(v_l_tlv[W_pair[1]]))
        Wlep2_dphi = LorentzVectorHEP.phi_mpi_pi(Z_phi - phi(v_l_tlv[W_pair[2]]))
        pt_4l = pt(sum(v_l_tlv)) / 1000
        if SR == 0
            sr_SF_inZ = 1
            sr_SF_noZ = 0
            sr_DF = 0
        elseif SR == 1
            sr_SF_inZ = 0
            sr_SF_noZ = 1
            sr_DF = 0
        else
            sr_SF_inZ = 0
            sr_SF_noZ = 0
            sr_DF = 1
        end

        NN_input = [HT, MET, METPhi, METSig, Njet, Wlep1_dphi, Wlep1_eta,
                    Wlep1_phi, Wlep1_pt, Wlep2_dphi, Wlep2_eta, Wlep2_phi,
                    Wlep2_pt, Zcand_mass, Zlep1_dphi, Zlep1_eta, Zlep1_phi,
                    Zlep1_pt, Zlep2_dphi, Zlep2_eta, Zlep2_phi, Zlep2_pt,
                    leptonic_HT, mass_4l, other_mass, pt_4l, total_HT,
                    sr_SF_inZ, sr_SF_noZ, sr_DF]
        NN_score = NN_calc(model, scales, minimums, NN_input)
        if controlregion == :ZZ
            SR != 0 && continue
            NN_score > 0.1 && continue
        end
        if !arrow_making
            @fill_dict! dict wgt atomic_push! pt_1, pt_2, pt_3, pt_4, eta_1, eta_2, 
            eta_3, eta_4, mass_4l, Zcand_mass, other_mass, METSig, MET, HT, leptonic_HT, total_HT,SR, 
            Z_eta, Z_phi, Z_pt, Z_rapidity, total_events, NN_score
        else
            jet_pt_1 = Njet < 1 ? 0.f0 : pt(v_j_tlv[v_j_order[1]]) / 1000
            jet_pt_2 = Njet < 2 ? 0.f0 : pt(v_j_tlv[v_j_order[2]]) / 1000
            jet_pt_3 = Njet < 3 ? 0.f0 : pt(v_j_tlv[v_j_order[3]]) / 1000
            jet_pt_4 = Njet < 4 ? 0.f0 : pt(v_j_tlv[v_j_order[4]]) / 1000
            jet_eta_1 = Njet < 1 ? -10.f0 : eta(v_j_tlv[v_j_order[1]])
            jet_eta_2 = Njet < 2 ? -10.f0 : eta(v_j_tlv[v_j_order[2]])
            jet_eta_3 = Njet < 3 ? -10.f0 : eta(v_j_tlv[v_j_order[3]])
            jet_eta_4 = Njet < 4 ? -10.f0 : eta(v_j_tlv[v_j_order[4]])
            jet_phi_1 = Njet < 1 ? -10.f0 : phi(v_j_tlv[v_j_order[1]])
            jet_phi_2 = Njet < 2 ? -10.f0 : phi(v_j_tlv[v_j_order[2]])
            jet_phi_3 = Njet < 3 ? -10.f0 : phi(v_j_tlv[v_j_order[3]])
            jet_phi_4 = Njet < 4 ? -10.f0 : phi(v_j_tlv[v_j_order[4]])
            jet_m_1 = Njet < 1 ? 0.f0 : mass(v_j_tlv[v_j_order[1]]) / 1000
            jet_m_2 = Njet < 2 ? 0.f0 : mass(v_j_tlv[v_j_order[2]]) / 1000
            jet_m_3 = Njet < 3 ? 0.f0 : mass(v_j_tlv[v_j_order[3]]) / 1000
            jet_m_4 = Njet < 4 ? 0.f0 : mass(v_j_tlv[v_j_order[4]]) / 1000
            (; v_j_btag60, v_j_btag70, v_j_btag77, v_j_btag85, v_j_btagCont) = evt
            jet_btagCont_1 = Njet < 1 ? -2 : v_j_btagCont[1]
            jet_btagCont_2 = Njet < 2 ? -2 : v_j_btagCont[2]
            jet_btagCont_3 = Njet < 3 ? -2 : v_j_btagCont[3]
            jet_btagCont_4 = Njet < 4 ? -2 : v_j_btagCont[4]
            @fill_dict! dict push! SR, Nlep, lep1_pid, lep2_pid, lep3_pid, lep4_pid, pt_1, pt_2, pt_3, pt_4, eta_1,
            eta_2, eta_3, eta_4, phi_1, phi_2, phi_3, phi_4, Njet, mass_4l, Zcand_mass, other_mass, MET,
            METSig, METPhi, leptonic_HT, HT, total_HT, Zlep1_pt, Zlep1_eta, Zlep1_phi, Zlep1_dphi, Zlep1_pid,
            Zlep2_pt, Zlep2_eta, Zlep2_phi, Zlep2_dphi, Zlep2_pid, Wlep1_pt, Wlep1_eta, Wlep1_phi, Wlep1_dphi,
            Wlep1_pid, Wlep2_pt, Wlep2_eta, Wlep2_phi, Wlep2_dphi, Wlep2_pid, chisq, pt_4l, jet_pt_1,
            jet_pt_2, jet_pt_3, jet_pt_4, jet_eta_1, jet_eta_2, jet_eta_3, jet_eta_4, jet_phi_1, jet_phi_2,
            jet_phi_3, jet_phi_4, jet_m_1, jet_m_2, jet_m_3, jet_m_4, v_j_btagCont, v_j_btag60,
            v_j_btag70, v_j_btag77, v_j_btag85, jet_btagCont_1, jet_btagCont_2, jet_btagCont_3, jet_btagCont_4, wgt
        end
    end
    return dict
end
