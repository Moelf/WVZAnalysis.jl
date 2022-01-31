function main_looper(mytree, sumWeight)
    # hists_dict = dictionary([
    #     :Z_mass_first => Hist1D(Float32; bins=0:10:200),
    #     :WZZ_ZZ_mass => Hist1D(Float32; bins=0:10:800),
    #     :WWZ_MET => Hist1D(Float32; bins=0:5:400),
    # ])

    data_ML = Dictionary(Dict(
        :SR => Int32[],
        :Nlep => Int32[],
#         :Zlep1_iso => Vector{Bool}[],
#         :Zlep2_iso => Vector{Bool}[],
#         :Wlep1_iso => Vector{Bool}[],
#         :Wlep2_iso => Vector{Bool}[],
#         :Zlep1_tight => Bool[],
#         :Zlep2_tight => Bool[],
#         :Wlep1_tight => Bool[],
#         :Wlep2_tight => Bool[],
#         :Zlep1_medium => Bool[],
#         :Zlep2_medium => Bool[],
#         :Wlep1_medium => Bool[],
#         :Wlep2_medium => Bool[],
        :pt_1 => Float32[],
        :pt_2 => Float32[],
        :pt_3 => Float32[],
        :pt_4 => Float32[],
        :Njet => Int32[],
        :mass_4l => Float32[],
        :Zcand_mass => Float32[],
        :other_mass => Float32[],
        :MET => Float32[],
        :METSig => Float32[],
        :METPhi => Float32[],
        :leptonic_HT => Float32[],
        :HT => Float32[],
        :total_HT => Float32[],
        :Zlep1_pt => Float32[],
        :Zlep1_eta => Float32[],
        :Zlep1_phi => Float32[],
        :Zlep1_dphi => Float32[],
        :Zlep1_mt => Float32[],
        :Zlep1_pid => Int32[],
        :Zlep2_pt => Float32[],
        :Zlep2_eta => Float32[],
        :Zlep2_phi => Float32[],
        :Zlep2_dphi => Float32[],
        :Zlep2_mt => Float32[],
        :Zlep2_pid => Int32[],
        :Wlep1_pt => Float32[],
        :Wlep1_eta => Float32[],
        :Wlep1_phi => Float32[],
        :Wlep1_dphi => Float32[],
        :Wlep1_mt => Float32[],
        :Wlep1_pid => Int32[],
        :Wlep2_pt => Float32[],
        :Wlep2_eta => Float32[],
        :Wlep2_phi => Float32[],
        :Wlep2_dphi => Float32[],
        :Wlep2_mt => Float32[],
        :Wlep2_pid => Int32[],
        :chisq => Float32[],
        :pt_4l => Float32[],
        :wgt => Float64[],
       ))
    @inbounds for (i, evt) in enumerate(mytree)
        
#         if (i-1 != 4)
#             continue
#         end
        
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
        other = setdiff((1,2,3,4), zpr1)
        wgt = evt.weight / sumWeight
        # push!(
        #     hists_dict[:Z_mass_first], mass(v_l_tlv[zpr1[1]] + v_l_tlv[zpr1[2]]) / 1000, wgt
        # )

        best_Z_mass = mass(v_l_tlv[zpr1[1]] + v_l_tlv[zpr1[2]])
        Z_phi = phi(v_l_tlv[zpr1[1]] + v_l_tlv[zpr1[2]])
        other_pair_mass = mass(v_l_tlv[other[1]] + v_l_tlv[other[2]])

        abs(best_Z_mass - Z_m) > 20e3 && continue

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
            continue
        end

        v_l_tight = vcat(evt.v_e_LHTight[e_mask], evt.v_m_tight[m_mask])
        v_l_medium = vcat(evt.v_e_LHMedium[e_mask], evt.v_m_medium[m_mask])
        pass_WWZ_cut, wgt, chi2, W_id = WWZ_Cut(
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
            l1, l2 = zpr1
            l3, l4 = W_id
            v_j_tlv = evt.v_j_tlv
            Njet = length(v_j_tlv)
            # Here we calculate the HT's:
            HT = isempty(v_j_tlv) ? 0.f0 : sum(pt, v_j_tlv) # hadronic HT
            leptonic_HT = sum([pt(v_l_tlv[v_l_order[x]]) for x in 1:4])
            total_HT    = HT + leptonic_HT


            # Putting a label to identify signal region:
            if (abs(v_l_pid[other[1]]) == abs(v_l_pid[other[2]]))
                if abs(other_pair_mass - Z_m) < 20e3
                    push!(data_ML[:SR], 0) # SF_inZ
                else
                    push!(data_ML[:SR], 1) # SF_noZ
                end
            else
                push!(data_ML[:SR], 2) # DF
            end

            push!(data_ML[:Nlep], length(v_l_pid))

#             push!(data_ML[:Zlep1_iso], v_l_passIso[l1])
#             push!(data_ML[:Zlep2_iso], v_l_passIso[l2])
#             push!(data_ML[:Wlep1_iso], v_l_passIso[l3])
#             push!(data_ML[:Wlep2_iso], v_l_passIso[l4])
#             push!(data_ML[:Zlep1_tight], v_l_tight[l1])
#             push!(data_ML[:Zlep2_tight], v_l_tight[l2])
#             push!(data_ML[:Wlep1_tight], v_l_tight[l3])
#             push!(data_ML[:Wlep2_tight], v_l_tight[l4])
#             push!(data_ML[:Zlep1_medium], v_l_medium[l1])
#             push!(data_ML[:Zlep2_medium], v_l_medium[l2])
#             push!(data_ML[:Wlep1_medium], v_l_medium[l3])
#             push!(data_ML[:Wlep2_medium], v_l_medium[l4])

            push!(data_ML[:Njet], Njet)
            push!(data_ML[:mass_4l], mass_4l)
            push!(data_ML[:Zcand_mass], best_Z_mass)
            push!(data_ML[:other_mass], other_pair_mass)
            push!(data_ML[:MET], evt.MET)
            push!(data_ML[:METSig], evt.METSig)
            push!(data_ML[:METPhi], evt.METPhi)
            push!(data_ML[:leptonic_HT], leptonic_HT)
            push!(data_ML[:HT], HT)
            push!(data_ML[:total_HT], total_HT)
            push!(data_ML[:pt_1], pt(v_l_tlv[v_l_order[1]]))
            push!(data_ML[:pt_2], pt(v_l_tlv[v_l_order[2]]))
            push!(data_ML[:pt_3], pt(v_l_tlv[v_l_order[3]]))
            push!(data_ML[:pt_4], pt(v_l_tlv[v_l_order[4]]))
            push!(data_ML[:Zlep1_pt], pt(v_l_tlv[l1]))
            push!(data_ML[:Zlep1_eta], eta(v_l_tlv[l1]))
            push!(data_ML[:Zlep1_phi], phi(v_l_tlv[l1]))
            push!(data_ML[:Zlep1_dphi], phi_mpi_pi(Z_phi - phi(v_l_tlv[l1])))
            push!(data_ML[:Zlep1_mt], mt(v_l_tlv[l1]))
            push!(data_ML[:Zlep1_pid], v_l_pid[l1])
            push!(data_ML[:Zlep2_pt], pt(v_l_tlv[l2]))
            push!(data_ML[:Zlep2_eta], eta(v_l_tlv[l2]))
            push!(data_ML[:Zlep2_phi], phi(v_l_tlv[l2]))
            push!(data_ML[:Zlep2_dphi], phi_mpi_pi(Z_phi - phi(v_l_tlv[l2])))
            push!(data_ML[:Zlep2_mt], mt(v_l_tlv[l2]))
            push!(data_ML[:Zlep2_pid], v_l_pid[l2])
            push!(data_ML[:Wlep1_pt], pt(v_l_tlv[l3]))
            push!(data_ML[:Wlep1_eta], eta(v_l_tlv[l3]))
            push!(data_ML[:Wlep1_phi], phi(v_l_tlv[l3]))
            push!(data_ML[:Wlep1_dphi], phi_mpi_pi(Z_phi - phi(v_l_tlv[l3])))
            push!(data_ML[:Wlep1_mt], mt(v_l_tlv[l3]))
            push!(data_ML[:Wlep1_pid], v_l_pid[l3])
            push!(data_ML[:Wlep2_pt], pt(v_l_tlv[l4]))
            push!(data_ML[:Wlep2_eta], eta(v_l_tlv[l4]))
            push!(data_ML[:Wlep2_phi], phi(v_l_tlv[l4]))
            push!(data_ML[:Wlep2_dphi], phi_mpi_pi(Z_phi - phi(v_l_tlv[l4])))
            push!(data_ML[:Wlep2_mt], mt(v_l_tlv[l4]))
            push!(data_ML[:Wlep2_pid], v_l_pid[l4])
            push!(data_ML[:chisq], chi2)
            push!(data_ML[:pt_4l], pt(sum(@view v_l_tlv[[l1,l2,l3,l4]])))
            # TODO HT
            push!(data_ML[:wgt], wgt)

            # push!(hists_dict[:WWZ_MET], evt.MET / 1000, wgt)
            continue
        end
    end
    # return hists_dict
    return data_ML
end