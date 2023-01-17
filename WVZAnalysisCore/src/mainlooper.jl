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
     shape_variation, sfsys) = task

    dict, pusher! = if arrow_making
        arrow_init(), push!
    elseif NN_hist
        NN_hist_init(; sfsys, shape_variation), push!
    else
        kinematic_hist_init(), push!
    end
    mytree = LazyTree(path, "tree_" * shape_variation)

    models = arrow_making ? nothing : init_BDT()
    main_looper(mytree, sumWeight, dict, models,
                shape_variation, sfsys, NN_hist, arrow_making, isdata)
end
function cutflow_total!(cutflow_ptr, dict, wgt_dict; NN_hist, shape_variation)
    if NN_hist && shape_variation == "NOMINAL"
        cutflow_ptr[] += 1
        push!(dict[:CutFlow], cutflow_ptr[])
        push!(dict[:CutFlowWgt], cutflow_ptr[], wgt_dict[:NOMINAL])
    end
end
function cutflow_SRs!(cutflow_ptr, dict, wgt_dict; NN_hist, SR, shape_variation)
    if NN_hist && shape_variation == "NOMINAL"
        cutflow_ptr[] += 1
        if SR == 0
            push!(dict[:SFinZCutFlow], cutflow_ptr[])
            push!(dict[:SFinZCutFlowWgt], cutflow_ptr[], wgt_dict[:NOMINAL])
        elseif SR == 1
            push!(dict[:SFnoZCutFlow], cutflow_ptr[])
            push!(dict[:SFnoZCutFlowWgt], cutflow_ptr[], wgt_dict[:NOMINAL])
        else
            push!(dict[:DFCutFlow], cutflow_ptr[])
            push!(dict[:DFCutFlowWgt], cutflow_ptr[], wgt_dict[:NOMINAL])
        end
    end
end

# above is function barrier
function main_looper(mytree, sumWeight, dict, models, 
        shape_variation, sfsys, NN_hist, arrow_making, isdata)
    model = models # for BDT
    for evt in mytree
        wgt_dict = Dict(:NOMINAL => 1 / sumWeight)
        if isdata
            wgt_dict[:NOMINAL] = 1.0
        else
            make_sfsys_wgt!(evt, wgt_dict,
                            :weight; sfsys, pre_mask=1)
            make_sfsys_wgt!(evt, wgt_dict, 
                            :v_j_wgt_btag77, Colon() ; sfsys)
            make_sfsys_wgt!(evt, wgt_dict, 
                            :w_sf_fjvt; sfsys, pre_mask=1)
        end
        ### initial_cut
        cutflow_ptr = Ref(0)
        cutflow_total!(cutflow_ptr, dict, wgt_dict; NN_hist, shape_variation)

        !(evt.passTrig) && continue
        cutflow_total!(cutflow_ptr, dict, wgt_dict; NN_hist, shape_variation)

        v_m_eta_orig, v_e_caloeta_orig = evt.v_m_eta, evt.v_e_caloeta
        e_etamask = [abs(η) < 2.47 && (abs(η)<1.37 || abs(η)>1.52) for η in v_e_caloeta_orig]
        m_etamask = [abs(η) < 2.5 for η in v_m_eta_orig]
        v_l_pid = @views Vcat(evt.v_e_pid[e_etamask], evt.v_m_pid[m_etamask])
        Nlep = length(v_l_pid)
        Nlep != 4 && continue
        if !isdata
            make_sfsys_wgt!(evt, wgt_dict,
                            :v_e_wgtReco, e_etamask; sfsys)
            make_sfsys_wgt!(evt, wgt_dict,
                            :v_m_wgtTTVA, m_etamask; sfsys)
        end
        cutflow_total!(cutflow_ptr, dict, wgt_dict; NN_hist, shape_variation)

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

        
        l3, l4 = W_pair
        sr_SF_inZ, sr_SF_noZ, sr_DF = if abs(v_l_pid[l3]) != abs(v_l_pid[l4]) #DF
            false, false, true
        elseif abs(other_mass - Z_m) < 20 # SF_inZ, already GeV
            true, false, false
        else #SF_noZ
            false, true, false
        end
        SR = sr_SF_inZ ? 0 : (sr_SF_noZ ? 1 : 2)
        cutflow_SRs!(cutflow_ptr, dict, wgt_dict; NN_hist, SR, shape_variation)

        mass_4l = mass(sum(v_l_tlv))
        mass_4l < 0.0 && continue
        ### end of initial_cut
        
        v_l_order = sortperm(v_l_tlv; by=pt, rev=true)

        ### QUALITY
        v_l_medium = @views Vcat(evt.v_e_LHMedium[e_etamask] , evt.v_m_medium[m_etamask])
        # for W leptons, require medium quality
        !v_l_medium[W_pair[1]] && continue
        !v_l_medium[W_pair[2]] && continue
        # for Z leptons only Loose requirement, no additional quality cut
        #
        if !isdata
            # I hate indexing
            Z_pair_in_e = filter(<=(Nelec), Z_pair)
            Z_pair_in_m = filter!(>(0), Z_pair .- Nelec)
            W_pair_in_e = filter(<=(Nelec), W_pair)
            W_pair_in_m = filter!(>(0), W_pair .- Nelec)
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
        end
        cutflow_SRs!(cutflow_ptr, dict, wgt_dict; NN_hist, SR, shape_variation)

        ### ISOLATION
        # for W leptons, require PLIVTight Isolation
        v_l_PLTight = @views Vcat(evt.v_e_passIso_PLImprovedTight[e_etamask], evt.v_m_passIso_PLImprovedTight[m_etamask])
        !v_l_PLTight[W_pair[1]] && continue
        !v_l_PLTight[W_pair[2]] && continue
        # for Z leptons require Loose isolation
        (;v_e_passIso_Loose_VarRad, v_m_passIso_PflowLoose_VarRad) = evt
        v_l_passIso = @views Vcat(v_e_passIso_Loose_VarRad[e_etamask], v_m_passIso_PflowLoose_VarRad[m_etamask])
        !v_l_passIso[Z_pair[1]] && continue
        !v_l_passIso[Z_pair[2]] && continue

        if !isdata
            make_sfsys_wgt!(evt, wgt_dict,
                            :v_e_wgtIso_Loose_VarRad_LooseBLayer, 
                            Z_pair_in_e; pre_mask = e_etamask, sfsys)
            make_sfsys_wgt!(evt, wgt_dict,
                            :v_m_wgtIso_PflowLoose_VarRad,
                            Z_pair_in_m; pre_mask = m_etamask, sfsys)
            make_sfsys_wgt!(evt, wgt_dict,
                            :v_e_wgtIso_PLImprovedTight_Medium, 
                            W_pair_in_e; pre_mask = e_etamask, sfsys)
            make_sfsys_wgt!(evt, wgt_dict,
                            :v_m_wgtIso_PLImprovedTight, 
                            W_pair_in_m; pre_mask = m_etamask, sfsys)
        end
        cutflow_SRs!(cutflow_ptr, dict, wgt_dict; NN_hist, SR, shape_variation)

        pass_mll = true
        chisq = WWZ_chi2(Z_pair, W_pair, v_l_pid, v_l_tlv)
        @inbounds for i in eachindex(v_l_tlv), j in (i + 1):length(v_l_tlv)
            v_l_pid[i] + v_l_pid[j] != 0 && continue
            WWZ_dilepton_mass = mass(v_l_tlv[i] + v_l_tlv[j])
            WWZ_dilepton_mass < 12 && (pass_mll = false)
        end
        !pass_mll && continue
        cutflow_SRs!(cutflow_ptr, dict, wgt_dict; NN_hist, SR, shape_variation)

        @inbounds pt(v_l_tlv[v_l_order[1]]) < 30 && continue
        @inbounds pt(v_l_tlv[v_l_order[2]]) < 15 && continue
        @inbounds pt(v_l_tlv[v_l_order[3]]) < 8 &&  continue
        @inbounds pt(v_l_tlv[v_l_order[4]]) < 6 &&  continue
        cutflow_SRs!(cutflow_ptr, dict, wgt_dict; NN_hist, SR, shape_variation)

        # selected lepton min dR
        pass_dR = true
        @inbounds for i in 1:4, j in (i + 1):4
            dR = deltar(v_l_tlv[v_l_order[i]], v_l_tlv[v_l_order[j]])
            dR < 0.1 && (pass_dR = false)
        end
        !pass_dR && continue
        cutflow_SRs!(cutflow_ptr, dict, wgt_dict; NN_hist, SR, shape_variation)
        
        has_b = any(evt.v_j_btag77)

        (; MET, METSig, METPhi) = evt
        MET /= 1000

        # force wgt to 1 for data


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
        Zlep1_dphi = LorentzVectorHEP.phi_mpi_pi(Z_phi - Zlep1_phi)
        Zlep2_dphi = LorentzVectorHEP.phi_mpi_pi(Z_phi - Zlep2_phi)
        Wlep1_pt, Wlep2_pt = pt.(@view v_l_tlv[W_pair])
        Wlep1_eta, Wlep2_eta = eta.(@view v_l_tlv[W_pair])
        Wlep1_phi, Wlep2_phi = phi.(@view v_l_tlv[W_pair])
        Wlep1_pid, Wlep2_pid = @view v_l_pid[W_pair]
        Wlep1_dphi = LorentzVectorHEP.phi_mpi_pi(Z_phi - Wlep1_phi)
        Wlep2_dphi = LorentzVectorHEP.phi_mpi_pi(Z_phi - Wlep2_phi)
        Wleps_deta = abs(Wlep1_eta - Wlep2_eta)
        pt_4l = pt(sum(v_l_tlv))

        MET_dPhi = LorentzVectorHEP.phi_mpi_pi(Z_phi - METPhi)
        
        event = evt.event
        # WARNING: don't use `mod1` it's shifting in the opposite direction
        moded_event = mod(event, 5) + 1
        if !arrow_making && NN_hist
            BDT_input = Float32[leptonic_HT, MET, Zlep1_dphi,
              Wlep1_pt, total_HT, Zlep2_dphi, Zlep2_eta, Njet, Wlep2_eta,
              Zlep2_pt, METSig, other_mass, Wlep1_dphi, Zlep1_pt,
              mass_4l, pt_4l, Zlep1_eta, HT, Wlep1_eta,
              Wlep2_dphi, Zcand_mass, Wlep2_pt, MET_dPhi]

            # NN_score = NN_calc(model, scales, minimums, NN_input)
            region = sr_SF_inZ ? :SFinZ : sr_SF_noZ ? :SFnoZ : :DF
            NN_score = model(BDT_input; fold = moded_event, region)
        end

        cr_ZZ = sr_SF_inZ && MET < 10 && !has_b
        cr_ttZ = has_b
        if cr_ZZ || cr_ttZ
            SR = -1
        end
        cutflow_ptr[] += 1

        if MET > 10 && NN_hist && shape_variation == "NOMINAL"
            if sr_SF_inZ push!(dict[:SFinZCutFlow], cutflow_ptr[]) end
            if sr_SF_inZ push!(dict[:SFinZCutFlowWgt], cutflow_ptr[], wgt_dict[:NOMINAL]) end
            if sr_SF_noZ push!(dict[:SFnoZCutFlow], cutflow_ptr[]) end
            if sr_SF_noZ push!(dict[:SFnoZCutFlowWgt], cutflow_ptr[], wgt_dict[:NOMINAL]) end
            if sr_DF push!(dict[:DFCutFlow], cutflow_ptr[]) end
            if sr_DF push!(dict[:DFCutFlowWgt], cutflow_ptr[], wgt_dict[:NOMINAL]) end
        end

        if NN_hist && !arrow_making
            region_prefix = if sr_SF_inZ
                :SFinZ
            elseif sr_SF_noZ
                :SFnoZ
            else
                :DF
            end
            if SR >= 0 && shape_variation == "NOMINAL"
                cutflow_ptr[] += 1
                if sr_SF_inZ push!(dict[:SFinZCutFlow], cutflow_ptr[]) end
                if sr_SF_inZ push!(dict[:SFinZCutFlowWgt], cutflow_ptr[], wgt_dict[:NOMINAL]) end
                if sr_SF_noZ push!(dict[:SFnoZCutFlow], cutflow_ptr[]) end
                if sr_SF_noZ push!(dict[:SFnoZCutFlowWgt], cutflow_ptr[], wgt_dict[:NOMINAL]) end
                if sr_DF push!(dict[:DFCutFlow], cutflow_ptr[]) end
                if sr_DF push!(dict[:DFCutFlowWgt], cutflow_ptr[], wgt_dict[:NOMINAL]) end
            end
            for (k,v) in wgt_dict
                if k==:NOMINAL && shape_variation != "NOMINAL"
                    k = shape_variation
                end
                if cr_ZZ
                    push!(dict[Symbol(:ZZCR, :__Njet, :__, k)], Njet, v)
                elseif cr_ttZ
                    push!(dict[Symbol(:ttZCR, :__Njet, :__, k)], Njet, v)
                elseif SR >= 0
                    push!(dict[Symbol(region_prefix, :__BDT, :__, k)], NN_score, v)
                    push!(dict[Symbol(region_prefix, :__MET, :__, k)], MET, v)
                    push!(dict[Symbol(region_prefix, :__Njet, :__, k)], Njet, v)
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
            @fill_dict! dict push! SR, Nlep, lep1_pid, lep2_pid, lep3_pid, lep4_pid, pt_1, pt_2, pt_3, pt_4, eta_1,
            eta_2, eta_3, eta_4, phi_1, phi_2, phi_3, phi_4, Njet, mass_4l, Zcand_mass, other_mass, MET,
            METSig, METPhi, MET_dPhi, leptonic_HT, HT, total_HT, Zlep1_pt, Zlep1_eta, Zlep1_phi, Zlep1_dphi, Zlep1_pid,
            Zlep2_pt, Zlep2_eta, Zlep2_phi, Zlep2_dphi, Zlep2_pid, Wlep1_pt, Wlep1_eta, Wlep1_phi, Wlep1_dphi,
            Wlep1_pid, Wlep2_pt, Wlep2_eta, Wlep2_phi, Wlep2_dphi, Wlep2_pid, Wleps_deta, chisq, pt_4l, jet_pt_1,
            jet_pt_2, jet_pt_3, jet_pt_4, jet_eta_1, jet_eta_2, jet_eta_3, jet_eta_4, jet_phi_1, jet_phi_2,
            jet_phi_3, jet_phi_4, jet_m_1, jet_m_2, jet_m_3, jet_m_4, v_j_btagCont, v_j_btag60,
            v_j_btag70, v_j_btag77, v_j_btag85, jet_btagCont_1, jet_btagCont_2, jet_btagCont_3, jet_btagCont_4, wgt, mcGenWgt,
            sr_SF_inZ, sr_SF_noZ, sr_DF, cr_ZZ, cr_ttZ, event
        else
            wgt = wgt_dict[:NOMINAL]
            @fill_dict! dict wgt push! pt_1, pt_2, pt_3, pt_4, eta_1, eta_2, 
            eta_3, eta_4, mass_4l, Zcand_mass, other_mass, METSig, MET, HT, leptonic_HT, total_HT, SR,
            Z_eta, Z_phi, Z_pt, Z_rapidity, Njet
        end
    end
    return dict
end