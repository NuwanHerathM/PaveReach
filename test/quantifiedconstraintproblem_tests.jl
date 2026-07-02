include("../src/quantifiedconstraintproblem.jl")

@testset "Get indices" begin
    sizes = [2, 3, 4]
    @test collect(build_ranges(sizes)) == [1:2, 3:5, 6:9]
    M = parseformula("P ∧ ¬ Q ∨ R")
    @test build_ranges_dict(M, sizes) == Dict("P" => 1:2, "Q" => 3:5, "R" => 6:9)
    @test build_positions_dict(M) == Dict("P" => 1, "Q" => 2, "R" => 3)
end