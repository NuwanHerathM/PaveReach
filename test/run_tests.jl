using Test

@testset "My tests" begin

    @testset "Pave tests" begin
        include("pave_tests.jl")
    end

    @testset "QuantifiedConstraintProblem tests" begin
        include("quantifiedconstraintproblem_tests.jl")
    end

    @testset "Utils tests" begin
        include("utils_tests.jl")
    end
end