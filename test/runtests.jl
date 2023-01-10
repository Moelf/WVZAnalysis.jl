using WVZAnalysis
using WVZAnalysis.FHist
using Test

if endswith(gethostname(), "uchicago.edu")
    @testset "WVZAnalysis.jl" begin
        fs = prep_tasks("Signal"; shape_variation="NOMINAL", sfsys=false);
        @test length(fs) == 33
        res = main_looper(sort(fs; by = x->filesize(x.path))[8]);
        @test bincounts(res[:SFnoZCutFlow]) == [0, 0, 0, 35, 19, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        @test integral(res[:SFnoZ__BDT__NOMINAL]) ≈ 0.0005643258391566838
    end
end
