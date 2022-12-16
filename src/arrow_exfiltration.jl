function arrow_init()
    return Dict(
                           #lep 1,2,3,4 ordered by decending Pt 
                           :SR => Int32[],
                           :event => UInt64[],
                           :sr_SF_inZ => Int32[],
                           :sr_SF_noZ => Int32[],
                           :sr_DF => Int32[],
                           :cr_ZZ => Int32[],
                           :cr_ttZ => Int32[],
                           :Nlep => Int32[],
                           :lep1_pid => Int32[], 
                           :lep2_pid => Int32[],
                           :lep3_pid => Int32[],
                           :lep4_pid => Int32[],
                           :pt_1 => Float32[],
                           :pt_2 => Float32[],
                           :pt_3 => Float32[],
                           :pt_4 => Float32[],
                           :eta_1 => Float32[],
                           :eta_2 => Float32[],
                           :eta_3 => Float32[],
                           :eta_4 => Float32[],
                           :phi_1 => Float32[],
                           :phi_2 => Float32[],
                           :phi_3 => Float32[],
                           :phi_4 => Float32[],
                           :Njet => Int32[],
                           :mass_4l => Float32[],
                           :Zcand_mass => Float32[],
                           :other_mass => Float32[],
                           :MET => Float32[],
                           :METSig => Float32[],
                           :METPhi => Float32[],
                           :MET_dPhi => Float32[],
                           :leptonic_HT => Float32[],
                           :HT => Float32[],
                           :total_HT => Float32[],
                           :Zlep1_pt => Float32[],
                           :Zlep1_eta => Float32[],
                           :Zlep1_phi => Float32[],
                           :Zlep1_dphi => Float32[],
                           :Zlep1_pid => Int32[],
                           :Zlep2_pt => Float32[],
                           :Zlep2_eta => Float32[],
                           :Zlep2_phi => Float32[],
                           :Zlep2_dphi => Float32[],
                           :Zlep2_pid => Int32[],
                           :Wlep1_pt => Float32[],
                           :Wlep1_eta => Float32[],
                           :Wlep1_phi => Float32[],
                           :Wlep1_dphi => Float32[],
                           :Wlep1_pid => Int32[],
                           :Wlep2_pt => Float32[],
                           :Wlep2_eta => Float32[],
                           :Wlep2_phi => Float32[],
                           :Wlep2_dphi => Float32[],
                           :Wleps_deta => Float32[],
                           :Wlep2_pid => Int32[],
                           :chisq => Float32[],
                           :pt_4l => Float32[],
                           :wgt => Float64[],
                           :mcGenWgt => Float32[],
                           :jet_pt_1 => Float32[],
                           :jet_pt_2 => Float32[],
                           :jet_pt_3 => Float32[],
                           :jet_pt_4 => Float32[],
                           :jet_eta_1 => Float32[],
                           :jet_eta_2 => Float32[],
                           :jet_eta_3 => Float32[],
                           :jet_eta_4 => Float32[],
                           :jet_phi_1 => Float32[],
                           :jet_phi_2 => Float32[],
                           :jet_phi_3 => Float32[],
                           :jet_phi_4 => Float32[],
                           :jet_m_1 => Float32[],
                           :jet_m_2 => Float32[],
                           :jet_m_3 => Float32[],
                           :jet_m_4 => Float32[],
                           :jet_m_4 => Float32[],
                           :v_j_btagCont => Vector{Int32}[],
                           :v_j_btag60 => Vector{Bool}[],
                           :v_j_btag70 => Vector{Bool}[],
                           :v_j_btag77 => Vector{Bool}[],
                           :v_j_btag85 => Vector{Bool}[],
                           :jet_btagCont_1 => Int32[],
                           :jet_btagCont_2 => Int32[],
                           :jet_btagCont_3 => Int32[],
                           :jet_btagCont_4 => Int32[],
                          )
end
