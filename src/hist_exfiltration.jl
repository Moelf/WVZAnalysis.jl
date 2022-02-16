@pour hist_prologue begin
  hists_dict = dictionary([
                           :Z_mass_first => Hist1D(Float32; bins=0:10:200),
                           :WZZ_ZZ_mass => Hist1D(Float32; bins=0:10:800),
                           :WWZ_MET => Hist1D(Float32; bins=0:5:400),
                          ])

end


@pour hist_epilogue begin

end
