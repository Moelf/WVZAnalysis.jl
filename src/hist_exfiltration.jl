@pour hist_prologue begin
  hists_dict = dictionary([
                           :WWZ_MET => Hist1D(Float32; bins=0:5:400),
                          ])

end


@pour hist_epilogue begin
    push!(hists_dict[:WWZ_MET], evt.MET)
end
