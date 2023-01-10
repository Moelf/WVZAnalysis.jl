using WVZAnalysisCore
using WVZAnalysisCore.FHist
using Test

@testset "WVZAnalysisCore.jl local test" begin
    path = joinpath(@__DIR__, "user.jiling.29896100._000001.ANALYSIS.root")
    t = AnalysisTask(;path, sumWeight=1)
    res = main_looper(t)
    @test bincounts(res[:CutFlow])[1:3] == [92, 84, 1]
    @test bincounts(res[:CutFlowWgt])[1:3] ≈ [18.961051007116527, 17.310845300169458, 0.23909674207864384]
end
# if endswith(gethostname(), "uchicago.edu")
#     @testset "WVZAnalysisCore.jl" begin
#         fs = prep_tasks("Signal"; shape_variation="NOMINAL", sfsys=false);
#         @test length(fs) == 33
#         res = main_looper(sort(fs; by = x->filesize(x.path))[8]);
#         @test bincounts(res[:SFnoZCutFlow]) == [0, 0, 0, 35, 19, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
#         @test integral(res[:SFnoZ__BDT__NOMINAL]) ≈ 0.0005643258391566838
#     end
# end
