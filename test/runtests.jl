using WVZAnalysis
using ThreadsX
using WVZAnalysis.FHist
using Test

@testset "single file test" begin
    path = joinpath(@__DIR__, "user.jiling.29896100._000001.ANALYSIS.root")
    t = AnalysisTask(;path, sumWeight=1)
    res = main_looper(t)
    @test bincounts(res[:CutFlow])[1:3] == [92, 84, 1]
    @test bincounts(res[:CutFlowWgt])[1:3] ≈ [19.637828743991196, 17.966735777895497, 0.23909674207864384]
end

@testset "prep_tasks dir test" begin
    tasks = prep_tasks("githubtest")
    t1 = main_looper(tasks[1])
    @test haskey(t1, :SFnoZ__BDT__NOMINAL)
end
