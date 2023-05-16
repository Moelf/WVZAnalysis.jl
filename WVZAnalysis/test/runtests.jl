using WVZAnalysis
using WVZAnalysis.FHist
using Test

@testset "WVZAnalysis.jl local test" begin
    path = joinpath(@__DIR__, "user.jiling.29896100._000001.ANALYSIS.root")
    t = AnalysisTask(;path, sumWeight=1)
    res = main_looper(t)
    @test bincounts(res[:CutFlow])[1:3] == [92, 84, 1]
    @test bincounts(res[:CutFlowWgt])[1:3] ≈ [19.637828743991196, 17.966735777895497, 0.23909674207864384]
end
