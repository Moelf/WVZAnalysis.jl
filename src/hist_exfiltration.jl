function kinematic_hist_init()
    return dictionary([
                           :pt_1 => Hist1D(Float64; bins=0:5:300, overflow=true),
                           :pt_2 => Hist1D(Float64; bins=0:5:300, overflow=true),
                           :pt_3 => Hist1D(Float64; bins=0:5:300, overflow=true),
                           :pt_4 => Hist1D(Float64; bins=0:5:300, overflow=true),
                           :eta_1 => Hist1D(Float64; bins=-3:0.2:3, overflow=true),
                           :eta_2 => Hist1D(Float64; bins=-3:0.2:3, overflow=true),
                           :eta_3 => Hist1D(Float64; bins=-3:0.2:3, overflow=true),
                           :eta_4 => Hist1D(Float64; bins=-3:0.2:3, overflow=true),
                           :mass_4l => Hist1D(Float64; bins=0:30:1000, overflow=true),
                           :Zcand_mass => Hist1D(Float64; bins=50:2:150, overflow=true),
                           :other_mass => Hist1D(Float64; bins=0:4:500, overflow=true),
                           :MET => Hist1D(Float64; bins=bins=0:1:1000, overflow=true),
                           :METSig => Hist1D(Float64; bins=0:0.2:20, overflow=true),
                           :HT => Hist1D(Float64; bins=0:10:800, overflow=true),
                           :leptonic_HT => Hist1D(Float64; bins=0:10:800, overflow=true),
                           :total_HT => Hist1D(Float64; bins=0:10:2000, overflow=true),
                           :SR => Hist1D(Float64; bins=0:3),
                           :Z_eta => Hist1D(Float64; bins=-4:0.4:4, overflow=true),
                           :Z_phi => Hist1D(Float64; bins=-4:0.4:4, overflow=true),
                           :Z_pt => Hist1D(Float64; bins=0:8:600, overflow=true),
                           :Z_rapidity => Hist1D(Float64; bins=-4:0.4:4, overflow=true),
                           :total_events => Hist1D(Float64; bins=0:1, overflow=true),
                           :NN_score => Hist1D(Float64; bins=0:0.01:1, overflow=true),
                          ])

end
