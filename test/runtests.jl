using WVZAnalysis
using Test

@testset "WVZAnalysis.jl" begin
    # Write your tests here.
    fs = prep_tasks("Signal"; shape_variation="NOMINAL", sfsys=false);
    @assert length(fs) == 33
    res = main_looper(fs[1])

end
