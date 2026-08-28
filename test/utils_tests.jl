include("../src/utils.jl")

@testset "ReLU derivative" begin
    @testset "Integers and floating point numbers" begin
        @test reluder(-1) == 0
        @test reluder(-0.256) == 0
        @test_throws DomainError reluder(0)
        @test reluder(1) == 1
        @test reluder(7.3) == 1
    end
    @testset "IntervalArithmetic" begin
        @test reluder(interval(-1, -0.5)) == interval(0, 0)
        @test reluder(interval(-1, 0)) == interval(0, 0)
        @test reluder(interval(-1, 1)) == interval(0, 1)
        @test reluder(interval(0, 1)) == interval(1, 1)
        @test reluder(interval(0.5, 1)) == interval(1, 1)
        @test reluder(interval(0, 0)) == interval(1, 1)
    end
end



@testset "Has a flat content" begin
    @test has_flat_content([1, 2, 3]) == true
    @test has_flat_content([1 2 3; 4 5 6]) == false
    @test has_flat_content([1]) == true
    @test has_flat_content([[[1]], [[2]]]) == true
    @test has_flat_content(reshape(1:6, (2, 3))) == false
    @test has_flat_content(reshape(1:6, (6, 1))) == true
    @test_throws AssertionError has_flat_content(Array{Int}(undef, 0, 3))
end