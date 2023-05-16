using WVZAnalysis
using WVZAnalysis.FHist
using Test

@testset "WVZAnalysis.jl local test" begin
    path = joinpath(@__DIR__, "user.jiling.29896100._000001.ANALYSIS.root")
    t = AnalysisTask(;path, sumWeight=1)
    res = main_looper(t)
    @test bincounts(res[:CutFlow])[1:3] == [92, 84, 1]
    @test bincounts(res[:CutFlowWgt])[1:3] ≈ [18.961051007116527, 17.310845300169458, 0.23909674207864384]
end
