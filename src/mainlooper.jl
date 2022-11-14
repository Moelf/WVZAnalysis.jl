"""
    main_looper(task::AnalysisTask)

The main entry point for running the main looper for a given task, it's done in steps:

1. destruct all the options from the `task` (see [`AnalysisTask`](@ref)):
```
    (; path, sumWeight, arrow_making, NN_hist, isdata, 
     shape_variation, controlregion, sfsys) = task
```
2. determine what output `dict` to prepare.
3. call the lower level `main_looper` (which has all arguments explicitly laied out)

This function also serves a [function barrier](https://docs.julialang.org/en/v1/manual/performance-tips/#kernel-functions)
for performance reason, because we have so many different behaviors in the same loop function.
"""
function main_looper(task::AnalysisTask)
    (; path, sumWeight, arrow_making, NN_hist, isdata, 
     shape_variation, controlregion, sfsys) = task

    dict, pusher! = if arrow_making
        arrow_init(), push!
    elseif NN_hist
        NN_hist_init(; sfsys, shape_variation), push!
    else
        kinematic_hist_init(), push!
    end
    mytree = LazyTree(path, "tree_" * shape_variation)

    models = arrow_making ? nothing : init_BDT()
    # models = init_ONNX()
    main_looper(mytree, sumWeight, dict, pusher!, models,
                shape_variation, sfsys, NN_hist, arrow_making, isdata, controlregion)
end

# above is function barrier
function main_looper(mytree, sumWeight, dict, pusher!, models, 
        shape_variation, sfsys, NN_hist, arrow_making, isdata, controlregion)
    model = models # for BDT
    # model, scales, minimums = models #for NN
    for evt in mytree
        ### initial_cut
        v_m_eta_orig, v_e_caloeta_orig = evt.v_m_eta, evt.v_e_caloeta
        e_etamask = [abs(η) < 2.47 && (abs(η)<1.37 || abs(η)>1.52) for η in v_e_caloeta_orig]
        m_etamask = [abs(η) < 2.5 for η in v_m_eta_orig]
        v_l_pid = @views Vcat(evt.v_e_pid[e_etamask], evt.v_m_pid[m_etamask])
        
        Nlep = length(v_l_pid)
        Nlep != 4 && continue
        Nelec = count(e_etamask)
        Nmuon = Nlep - Nelec
        (; v_e_pt, v_e_eta, v_e_phi,
         v_m_pt, v_m_eta, v_m_phi) = evt
        v_e_pt = v_e_pt ./ 1000
        v_m_pt = v_m_pt ./ 1000
        v_e_m = fill!(similar(v_e_pt), e_mass)
        v_m_m = fill!(similar(v_m_pt), m_mass)

        v_e_tlv = @views LorentzVectorCyl.(v_e_pt[e_etamask], v_e_eta[e_etamask], v_e_phi[e_etamask], v_e_m[e_etamask])
        v_m_tlv = @views LorentzVectorCyl.(v_m_pt[m_etamask], v_m_eta[m_etamask], v_m_phi[m_etamask], v_m_m[m_etamask])
        v_l_tlv = Vcat(v_e_tlv, v_m_tlv)
        Z_pair, W_pair, best_Z_mass = Find_Z_Pairs(v_l_pid, v_l_tlv)
        isinf(best_Z_mass) && continue
        other_mass = mass(v_l_tlv[W_pair[1]] + v_l_tlv[W_pair[2]])
        abs(best_Z_mass - Z_m) > 20 && continue
        mass_4l = mass(sum(v_l_tlv))
        mass_4l < 0.0 && continue
        ### end of initial_cut
        !(evt.passTrig) && continue
        v_l_order = sortperm(v_l_tlv; by=pt, rev=true)
        v_l_medium = @views Vcat(evt.v_e_LHMedium[e_etamask] , evt.v_m_medium[m_etamask]) #quality

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
        end
        (;
         v_e_passIso_Loose_VarRad,
         v_m_passIso_PflowLoose_VarRad,
        ) = evt
        v_l_passIso = @views Vcat(v_e_passIso_Loose_VarRad[e_etamask], v_m_passIso_PflowLoose_VarRad[m_etamask])

        pass_WWZ_cut, chisq, W_id = WWZ_Cut(
                                            Z_pair,
                                            W_pair,
                                            v_l_pid,
                                            v_l_order,
                                            v_l_tlv
                                           )
        for i in 1:2
            ### for Z leptons isolation: Loose(e) and Loose(mu)
            ( (abs(v_l_pid[Z_pair[i]]) == 11) && !v_l_passIso[Z_pair[i]] ) && (pass_WWZ_cut = false)
            ( (abs(v_l_pid[Z_pair[i]]) == 13) && !v_l_passIso[Z_pair[i]] ) && (pass_WWZ_cut = false)

            ### for W leptons, require medium quality and PLIV tight
            ( !v_l_medium[W_pair[i]] ) && (pass_WWZ_cut = false)
            isdata && break
        end

        !pass_WWZ_cut && continue
        
        b_idx, has_b = Bjet_cut(evt)
        (; MET, METSig, METPhi) = evt
        MET /= 1000
        if controlregion == :ttZ
            MET < 20 && continue
            (other_mass > 80 && other_mass < 100) && continue
            !has_b && continue
        else
            has_b && continue
        end
        wgt_dict = Dict(:NOMINAL => 1 / sumWeight)
        make_sfsys_wgt!(evt, wgt_dict,
                        :weight; sfsys, pre_mask=1)
        make_sfsys_wgt!(evt, wgt_dict, 
                        :v_j_wgt_btag77, b_idx; sfsys)

        l3, l4 = W_pair
        # force wgt to 1 for data
        if isdata
            wgt_dict[:NOMINAL] = 1.0
        else
            # I hate indexing
            Z_pair_in_e = filter(<=(Nelec), Z_pair)
            Z_pair_in_m = filter!(>(0), Z_pair .- Nelec)
            W_pair_in_e = filter(<=(Nelec), W_pair)
            W_pair_in_m = filter!(>(0), W_pair .- Nelec)

            make_sfsys_wgt!(evt, wgt_dict,
                            :v_e_wgtReco, e_etamask; sfsys)
            make_sfsys_wgt!(evt, wgt_dict,
                            :v_m_wgtTTVA, m_etamask; sfsys)
            if controlregion!=:ZJets
                # v_l_wgtPLTight = Vcat(evt.v_e_wgtIso_PLImprovedTight_Medium[e_etamask], 
                # evt.v_m_wgtIso_PLImprovedTight[m_etamask])
                # wgt *= @views reduce(*, v_l_wgtPLTight[W_pair]; init = 1.0)
                make_sfsys_wgt!(evt, wgt_dict,
                                :v_e_wgtIso_PLImprovedTight_Medium, 
                                W_pair_in_e; pre_mask = e_etamask, sfsys)
                make_sfsys_wgt!(evt, wgt_dict,
                                :v_m_wgtIso_PLImprovedTight, 
                                W_pair_in_m; pre_mask = m_etamask, sfsys)
            end
            # quality weights
            # v_l_wgtLoose =  @views Vcat(evt.v_e_wgtLoose[e_etamask] , evt.v_m_wgtLoose[m_etamask]) # quality wgt
            # v_l_wgtMedium = @views Vcat(evt.v_e_wgtMedium[e_etamask], evt.v_m_wgtMedium[m_etamask]) # quality wgt
            make_sfsys_wgt!(evt, wgt_dict,
                            :v_e_wgtLoose,
                            Z_pair_in_e; pre_mask = e_etamask, sfsys)
            make_sfsys_wgt!(evt, wgt_dict,
                            :v_m_wgtLoose,
                            Z_pair_in_m; pre_mask = m_etamask, sfsys)
            make_sfsys_wgt!(evt, wgt_dict,
                            :v_e_wgtMedium,
                            W_pair_in_e; pre_mask = e_etamask, sfsys)
            make_sfsys_wgt!(evt, wgt_dict,
                            :v_m_wgtMedium,
                            W_pair_in_m; pre_mask = m_etamask, sfsys)
            # iso weights
            # v_l_wgtIso =  @views Vcat(v_e_wgtIso_Loose_VarRad_LooseBLayer[e_etamask], v_m_wgtIso_PflowLoose_VarRad[m_etamask])
            make_sfsys_wgt!(evt, wgt_dict,
                            :v_e_wgtIso_Loose_VarRad_LooseBLayer, 
                            Z_pair_in_e; pre_mask = e_etamask, sfsys)
            make_sfsys_wgt!(evt, wgt_dict,
                            :v_m_wgtIso_PflowLoose_VarRad,
                            Z_pair_in_m; pre_mask = m_etamask, sfsys)
        end

        lep1_pid, lep2_pid, lep3_pid, lep4_pid = @view v_l_pid[v_l_order]
        pt_1, pt_2, pt_3, pt_4 = pt.(@view v_l_tlv[v_l_order])
        eta_1, eta_2, eta_3, eta_4 = eta.(@view v_l_tlv[v_l_order])
        phi_1, phi_2, phi_3, phi_4 = phi.(@view v_l_tlv[v_l_order])
        phi_1, phi_2, phi_3, phi_4 = phi.(@view v_l_tlv[v_l_order])
        v_j_tlv = LorentzVectorCyl.(evt.v_j_pt ./ 1000, evt.v_j_eta, evt.v_j_phi, evt.v_j_m ./ 1000)
        v_j_order = sortperm(v_j_tlv; by=pt, rev=true)
        Njet = length(v_j_tlv)
        Zcand_mass = best_Z_mass
        Z_eta = eta(v_l_tlv[Z_pair[1]]+v_l_tlv[Z_pair[2]])
        Z_phi = phi(v_l_tlv[Z_pair[1]]+v_l_tlv[Z_pair[2]])
        Z_pt = pt(v_l_tlv[Z_pair[1]]+v_l_tlv[Z_pair[2]])
        HT = sum(pt, v_j_tlv; init=0.f0) # hadronic HT
        leptonic_HT = sum(pt, v_l_tlv; init=0.f0)
        total_HT    = HT + leptonic_HT
  
        Z_tlv = v_l_tlv[Z_pair[1]]+v_l_tlv[Z_pair[2]]
        Z_rapidity = 0.5 * log((energy(Z_tlv)+pz(Z_tlv))/(energy(Z_tlv)-pz(Z_tlv)))
        Zlep1_pt, Zlep2_pt = pt.(@view v_l_tlv[Z_pair])
        Zlep1_eta, Zlep2_eta = eta.(@view v_l_tlv[Z_pair])
        Zlep1_phi, Zlep2_phi = phi.(@view v_l_tlv[Z_pair])
        Zlep1_pid, Zlep2_pid = @view v_l_pid[Z_pair]
        Zlep1_dphi = LorentzVectorHEP.phi_mpi_pi(Z_phi - phi(v_l_tlv[Z_pair[1]]))
        Zlep2_dphi = LorentzVectorHEP.phi_mpi_pi(Z_phi - phi(v_l_tlv[Z_pair[2]]))
        Wlep1_pt, Wlep2_pt = pt.(@view v_l_tlv[W_pair])
        Wlep1_eta, Wlep2_eta = eta.(@view v_l_tlv[W_pair])
        Wlep1_phi, Wlep2_phi = phi.(@view v_l_tlv[W_pair])
        Wlep1_pid, Wlep2_pid = @view v_l_pid[W_pair]
        Wlep1_dphi = LorentzVectorHEP.phi_mpi_pi(Z_phi - phi(v_l_tlv[W_pair[1]]))
        Wlep2_dphi = LorentzVectorHEP.phi_mpi_pi(Z_phi - phi(v_l_tlv[W_pair[2]]))
        pt_4l = pt(sum(v_l_tlv))

        MET_dPhi = LorentzVectorHEP.phi_mpi_pi(Z_phi - METPhi)

        sr_SF_inZ, sr_SF_noZ, sr_DF = if abs(v_l_pid[l3]) != abs(v_l_pid[l4]) #DF
            false, false, true
        elseif abs(other_mass - Z_m) < 20 # SF_inZ, already GeV
            true, false, false
        else #SF_noZ
            false, true, false
        end
        SR = sr_SF_inZ ? 0 : (sr_SF_noZ ? 1 : 2)
        
        wgt = wgt_dict[:NOMINAL]
        if !arrow_making && NN_hist
            # NN_input = Float32[HT, MET, METPhi, METSig, Njet, Wlep1_dphi, Wlep1_eta,
            #             Wlep1_phi, Wlep1_pt, Wlep2_dphi, Wlep2_eta, Wlep2_phi,
            #             Wlep2_pt, Zcand_mass, Zlep1_dphi, Zlep1_eta, Zlep1_phi,
            #             Zlep1_pt, Zlep2_dphi, Zlep2_eta, Zlep2_phi, Zlep2_pt,
            #             leptonic_HT, mass_4l, other_mass, pt_4l, total_HT,
            #             sr_SF_inZ, sr_SF_noZ, sr_DF]
            BDT_input = Float32[Zlep1_phi, leptonic_HT, Zlep2_phi, MET, Zlep1_dphi,
              Wlep1_pt, total_HT, Zlep2_dphi, Zlep2_eta, Njet, Wlep2_eta,
              Zlep2_pt, METSig, other_mass, Wlep1_dphi, Zlep1_pt, METPhi,
              mass_4l, pt_4l, Wlep2_phi, Zlep1_eta, HT, Wlep1_eta,
              Wlep2_dphi, Zcand_mass, Wlep2_pt, Wlep1_phi, sr_SF_inZ,
              sr_SF_noZ, sr_DF]

            # NN_score = NN_calc(model, scales, minimums, NN_input)
            NN_score = model(BDT_input)
        end

        is_CR = begin
            sr_SF_inZ &&
            MET < 10 &&
            Njet == 0
        end

        if NN_hist && !arrow_making
            region_prefix = if sr_SF_inZ
                :SFinZ
            elseif sr_SF_noZ
                :SFnoZ
            else
                :DF
            end
            for (k,v) in wgt_dict
                if is_CR
                    pusher!(dict[Symbol(:CR, :__yield, :__, k)], NN_score, v)
                else
                    if k==:NOMINAL && shape_variation != "NOMINAL"
                        k = shape_variation
                    end
                    pusher!(dict[Symbol(region_prefix, :__NN, :__, k)], NN_score, v)
                    pusher!(dict[Symbol(region_prefix, :__MET, :__, k)], MET, v)
                end
            end

        elseif arrow_making
            jet_pt_1 = Njet < 1 ? 0.f0 : pt(v_j_tlv[v_j_order[1]])
            jet_pt_2 = Njet < 2 ? 0.f0 : pt(v_j_tlv[v_j_order[2]])
            jet_pt_3 = Njet < 3 ? 0.f0 : pt(v_j_tlv[v_j_order[3]])
            jet_pt_4 = Njet < 4 ? 0.f0 : pt(v_j_tlv[v_j_order[4]])
            jet_eta_1 = Njet < 1 ? -10.f0 : eta(v_j_tlv[v_j_order[1]])
            jet_eta_2 = Njet < 2 ? -10.f0 : eta(v_j_tlv[v_j_order[2]])
            jet_eta_3 = Njet < 3 ? -10.f0 : eta(v_j_tlv[v_j_order[3]])
            jet_eta_4 = Njet < 4 ? -10.f0 : eta(v_j_tlv[v_j_order[4]])
            jet_phi_1 = Njet < 1 ? -10.f0 : phi(v_j_tlv[v_j_order[1]])
            jet_phi_2 = Njet < 2 ? -10.f0 : phi(v_j_tlv[v_j_order[2]])
            jet_phi_3 = Njet < 3 ? -10.f0 : phi(v_j_tlv[v_j_order[3]])
            jet_phi_4 = Njet < 4 ? -10.f0 : phi(v_j_tlv[v_j_order[4]])
            jet_m_1 = Njet < 1 ? 0.f0 : mass(v_j_tlv[v_j_order[1]])
            jet_m_2 = Njet < 2 ? 0.f0 : mass(v_j_tlv[v_j_order[2]])
            jet_m_3 = Njet < 3 ? 0.f0 : mass(v_j_tlv[v_j_order[3]])
            jet_m_4 = Njet < 4 ? 0.f0 : mass(v_j_tlv[v_j_order[4]])
            (; v_j_btag60, v_j_btag70, v_j_btag77, v_j_btag85, v_j_btagCont) = evt
            jet_btagCont_1 = Njet < 1 ? -2 : v_j_btagCont[1]
            jet_btagCont_2 = Njet < 2 ? -2 : v_j_btagCont[2]
            jet_btagCont_3 = Njet < 3 ? -2 : v_j_btagCont[3]
            jet_btagCont_4 = Njet < 4 ? -2 : v_j_btagCont[4]
            mcGenWgt = if isdata
                1.0
            else
                first(evt.v_mcGenWgt)
            end
            event = evt.event
            @fill_dict! dict pusher! SR, Nlep, lep1_pid, lep2_pid, lep3_pid, lep4_pid, pt_1, pt_2, pt_3, pt_4, eta_1,
            eta_2, eta_3, eta_4, phi_1, phi_2, phi_3, phi_4, Njet, mass_4l, Zcand_mass, other_mass, MET,
            METSig, METPhi, MET_dPhi, leptonic_HT, HT, total_HT, Zlep1_pt, Zlep1_eta, Zlep1_phi, Zlep1_dphi, Zlep1_pid,
            Zlep2_pt, Zlep2_eta, Zlep2_phi, Zlep2_dphi, Zlep2_pid, Wlep1_pt, Wlep1_eta, Wlep1_phi, Wlep1_dphi,
            Wlep1_pid, Wlep2_pt, Wlep2_eta, Wlep2_phi, Wlep2_dphi, Wlep2_pid, chisq, pt_4l, jet_pt_1,
            jet_pt_2, jet_pt_3, jet_pt_4, jet_eta_1, jet_eta_2, jet_eta_3, jet_eta_4, jet_phi_1, jet_phi_2,
            jet_phi_3, jet_phi_4, jet_m_1, jet_m_2, jet_m_3, jet_m_4, v_j_btagCont, v_j_btag60,
            v_j_btag70, v_j_btag77, v_j_btag85, jet_btagCont_1, jet_btagCont_2, jet_btagCont_3, jet_btagCont_4, wgt, mcGenWgt,
            sr_SF_inZ, sr_SF_noZ, sr_DF, event
        else
            @fill_dict! dict wgt pusher! pt_1, pt_2, pt_3, pt_4, eta_1, eta_2, 
            eta_3, eta_4, mass_4l, Zcand_mass, other_mass, METSig, MET, HT, leptonic_HT, total_HT,SR, 
            Z_eta, Z_phi, Z_pt, Z_rapidity, Njet
        end
    end
    return dict
end
