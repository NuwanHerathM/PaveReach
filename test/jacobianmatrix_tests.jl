@testset "Jacobian" begin
    MJ = MatJacobianMatrix([(x...) -> x[1] (x...) -> 2*x[2]; (x...) -> (x[2] + 3) (x...) -> (x[1] - 1)])
    RJ = RowJacobianMatrix([(x...) -> [x[1], 2*x[2]], (x...) -> [(x[2] + 3), (x[1] - 1)]])
    CJ = ColJacobianMatrix([(x...) -> [x[1], (x[2] + 3)], (x...) -> [2*x[2], (x[1] - 1)]])
    
    @testset "Matrix equivalence" begin
        result = [1 2; 4 0]
        for i in 1:2
            for j in 1:2
                @test MJ[i, j](1 ,1) == RJ[i, j](1 ,1) == CJ[i, j](1 ,1) == result[i, j]
            end
        end
    end
    
    @testset "Gradient" begin
        result = [1, 2]
        @test get_gradient(MJ, 1)(1, 1) == get_gradient(CJ, 1)(1, 1) == get_gradient(RJ, 1)(1, 1) == result
    end

    @testset "Matrix inflation" begin
        MJ_inflated = inflate(MJ, [1])
        CJ_inflated = inflate(CJ, [1])
        result = vcat(get_gradient(MJ, 1), get_gradient(MJ, 1), get_gradient(MJ, 2))
        for i in 1:3
            for j in 1:2
                @test MJ_inflated[i, j](1, 1) == CJ_inflated[i, j](1, 1)
            end
        end
    end

    @testset "Append constant matrix" begin
        Z = Matrix{Float64}(-LinearAlgebra.I, 2, 2)
        MJZ = append_constant(MJ, Z)
        RJZ = append_constant(RJ, Z)
        CJZ = append_constant(CJ, Z)
        result = [1 2 -1 0; 4 0 0 -1]
        for i in 1:2
            for j in 1:4
                @test MJZ[i,j](1, 1) == result[i, j]
            end
            @test get_gradient(RJZ, i)(1, 1) == result[i, :]
        end
        for j in 1:4
            @test CJZ.matrix[1, j](1, 1) == result[:, j]
        end
    end
end