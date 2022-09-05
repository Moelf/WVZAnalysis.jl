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

    model, rescaling_parameters = init_ONNX()

    @floop executor for evt in mytree
    # for evt in mytree
        ### initial_cut
        pt_1 = 5
        v_m_eta_orig, v_e_caloeta_orig = evt.v_m_eta, evt.v_e_caloeta
        e_etamask = [abs(η) < 2.47 && (abs(η)<1.37 || abs(η)>1.52) for η in v_e_caloeta_orig]
        m_etamask = [abs(η) < 2.5 for η in v_m_eta_orig]
        #e_etamask = [true for η in v_e_caloeta_orig]
        #m_etamask = [true for η in v_m_eta_orig]
        v_l_pid = @views Vcat(evt.v_e_pid[e_etamask], evt.v_m_pid[m_etamask])
        
        wgt = evt.weight / sumWeight * wgt_factor

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
        v_l_order = sortperm(v_l_tlv; by=pt, rev=true)
        v_l_medium = @views Vcat(evt.v_e_LHMedium[e_etamask] , evt.v_m_medium[m_etamask]) #quality
        if isdata
            v_l_wgtLoose = [1 for lepton in v_l_tlv]
            v_l_wgtMedium = [1 for lepton in v_l_tlv]
        else
            v_l_wgtLoose =  @views Vcat(evt.v_e_wgtLoose[e_etamask] , evt.v_m_wgtLoose[m_etamask]) # quality wgt
            v_l_wgtMedium = @views Vcat(evt.v_e_wgtMedium[e_etamask], evt.v_m_wgtMedium[m_etamask]) #quality wgt
            wgt *= @views reduce(*, evt.v_e_wgtReco[e_etamask])
            wgt *= @views reduce(*, evt.v_m_wgtTTVA[m_etamask])
        end
        b_wgt, b_veto = Bjet_Cut(evt)
        wgt *= b_wgt

        cut = 0.5

        Nlep = length(v_l_pid)

        
        if length(v_m_eta_orig)+length(v_e_caloeta_orig) == 4 && Nlep != 4
            @fill_dict! dict wgt atomic_push! cut, pt_1
            continue
        end
        #if sum(e_etamask) != length(v_e_caloeta_orig) || sum(m_etamask) != length(v_m_eta_orig)
        #    @fill_dict! dict wgt atomic_push! cut, pt_1
        #    continue
        #end
        cut+=1
        
        if Nlep != 4
            @fill_dict! dict wgt atomic_push! cut, pt_1
            continue
        end
        cut+=1
        
        # force wgt to 1 for data
        if isdata
            wgt = 1.0
        end
        
        if isinf(best_Z_mass)
            @fill_dict! dict wgt atomic_push! cut, pt_1
            continue
        end
        cut+=1
        
        other_mass = mass(v_l_tlv[W_pair[1]] + v_l_tlv[W_pair[2]])

        SR = if abs(v_l_pid[W_pair[1]]) != abs(v_l_pid[W_pair[2]]) #DF
            2
        elseif abs(other_mass - Z_m) < 20e3 # SF_inZ
            0
        else #SF_noZ
            1
        end

        if abs(best_Z_mass - Z_m) > 20e3
            @fill_dict! dict wgt atomic_push! cut, pt_1
            atomic_push!(dict[:SRcut],SR,cut,wgt)
            continue
        end
        cut+=1

        mass_4l = mass(sum(v_l_tlv))

        if mass_4l < 0.0
            @fill_dict! dict wgt atomic_push! cut, pt_1
            atomic_push!(dict[:SRcut],SR,cut,wgt)
            continue
        end
        cut+=1

        ### end of initial_cut
        if !(evt.passTrig)
            @fill_dict! dict wgt atomic_push! cut, pt_1
            atomic_push!(dict[:SRcut],SR,cut,wgt)
            continue
        end
        cut+=1

        
        
        ############## use PLIV for W lepton ISO #################
        v_l_PLTight = Vcat(evt.v_e_passIso_PLImprovedTight[e_etamask], evt.v_m_passIso_PLImprovedTight[m_etamask])
        failed_PLTight = any(==(false), v_l_PLTight[W_pair])

        if failed_PLTight
            @fill_dict! dict wgt atomic_push! cut, pt_1
            atomic_push!(dict[:SRcut],SR,cut,wgt)
            continue
        end
        cut+=1

        if !isdata
            v_l_wgtPLTight = Vcat(evt.v_e_wgtIso_PLImprovedTight_Medium[e_etamask], evt.v_m_wgtIso_PLImprovedTight[m_etamask])
            wgt *= @views reduce(*, v_l_wgtPLTight[W_pair]; init = 1.0)
        end
        (;
         v_e_passIso_Loose_VarRad,
         v_m_passIso_PflowLoose_VarRad,
         v_e_wgtIso_Loose_VarRad_LooseBLayer,
         v_m_wgtIso_PflowLoose_VarRad,
        ) = evt
        v_l_passIso = @views Vcat(v_e_passIso_Loose_VarRad[e_etamask], v_m_passIso_PflowLoose_VarRad[m_etamask])
        v_l_wgtIso =  @views Vcat(v_e_wgtIso_Loose_VarRad_LooseBLayer[e_etamask], v_m_wgtIso_PflowLoose_VarRad[m_etamask])
        

        chi2 = WWZ_chi2(Z_pair, W_pair, v_l_pid, v_l_tlv)

        mllcut = false
        @inbounds for i in eachindex(v_l_tlv)
            for j in (i + 1):length(v_l_tlv)
                v_l_pid[i] + v_l_pid[j] != 0 && continue
                WWZ_dilepton_mass = mass(v_l_tlv[i] + v_l_tlv[j]) / 1000
                if WWZ_dilepton_mass < 12
                    mllcut = true
                end
            end
        end
        if mllcut
            @fill_dict! dict wgt atomic_push! cut, pt_1
            atomic_push!(dict[:SRcut],SR,cut,wgt)
            continue
        end
        cut+=1
        
        ptcut = false
        if @inbounds pt(v_l_tlv[v_l_order[1]]) < 30e3 || @inbounds pt(v_l_tlv[v_l_order[2]]) < 15e3 || @inbounds pt(v_l_tlv[v_l_order[3]]) < 8e3 || @inbounds pt(v_l_tlv[v_l_order[4]]) < 6e3
            ptcut = true
        end
        if ptcut
            @fill_dict! dict wgt atomic_push! cut, pt_1
            atomic_push!(dict[:SRcut],SR,cut,wgt)
            continue
        end
        cut+=1
        drcut = false
        # selected lepton min dR
        @inbounds for i in 1:4
            for j in (i + 1):4
                dR = deltar(v_l_tlv[v_l_order[i]], v_l_tlv[v_l_order[j]])
                if dR < 0.1
                    drcut = true
                end
            end
        end
        if drcut
            @fill_dict! dict wgt atomic_push! cut, pt_1
            atomic_push!(dict[:SRcut],SR,cut,wgt)
            continue
        end
        cut+=1

        WWZ_wgt = wgt
        Zisocut = false
        Wqualcut = false
        for i in 1:2
            ### for Z leptons isolation: Loose(e) and Loose(mu)
            if ( (abs(v_l_pid[Z_pair[i]]) == 11) && !v_l_passIso[Z_pair[i]] ) || ( (abs(v_l_pid[Z_pair[i]]) == 13) && !v_l_passIso[Z_pair[i]] )
                Zisocut = true
            end

            ### for W leptons, require medium quality and PLIV tight
            if ( !v_l_medium[W_pair[i]] )
                Wqualcut = true
            end
            isdata && continue
            # quality weights
            WWZ_wgt *= v_l_wgtLoose[Z_pair[i]] * v_l_wgtMedium[W_pair[i]]
            # iso weights
            WWZ_wgt *= v_l_wgtIso[Z_pair[i]]
        end
    
        if Zisocut
            @fill_dict! dict wgt atomic_push! cut, pt_1
            atomic_push!(dict[:SRcut],SR,cut,wgt)
            continue
        end
        cut+=1

        if Wqualcut
            @fill_dict! dict wgt atomic_push! cut, pt_1
            atomic_push!(dict[:SRcut],SR,cut,wgt)
            continue
        end
        cut+=1

        wgt = WWZ_wgt
        # in `b_veto == true` means we've passed the criterial,
        # which means we didn't see a b-jet
        
        if !b_veto
            @fill_dict! dict wgt atomic_push! cut, pt_1
            atomic_push!(dict[:SRcut],SR,cut,wgt)
            continue
        end
        cut+=1

        @fill_dict! dict wgt atomic_push! cut, pt_1
        atomic_push!(dict[:SRcut],SR,cut,wgt)
        
        
    end
    return dict
end