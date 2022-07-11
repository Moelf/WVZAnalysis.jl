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
                           :Zcand_mass => Hist1D(Float64; bins=70:2:110, overflow=true),
                           :other_mass => Hist1D(Float64; bins=0:5:200, overflow=true),
                           :MET => Hist1D(Float64; bins=0:5:300, overflow=true),
                           :METSig => Hist1D(Float64; bins=0:0.2:20, overflow=true),
                           :HT => Hist1D(Float64; bins=0:10:800, overflow=true),
                           :leptonic_HT => Hist1D(Float64; bins=0:10:800, overflow=true),
                           :total_HT => Hist1D(Float64; bins=0:10:2000, overflow=true),
                           :SR => Hist1D(Float64; bins=0:3),
                          ])

end
