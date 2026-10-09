abstract type JacobianMatrix end

# MatJacobianMatrix

"""
    Jacobian storing the M * N functions ∂f_j/∂x_i for f: R^N → R^M.
"""
struct MatJacobianMatrix{S1, S2} <: JacobianMatrix
    matrix::SMatrix{S1, S2, Function}
    function MatJacobianMatrix(matrix::AbstractMatrix)
        S1, S2 = size(matrix)
        new{S1, S2}(SMatrix{S1, S2, Function}(matrix))
    end
end

# RowJacobianMatrix

"""
    Jacobian storing the M functions ∂f_j/∂x for f: R^N → R^M, that is the gradients ∇f_j.
"""
struct RowJacobianMatrix{S} <: JacobianMatrix
    matrix::SMatrix{S, 1, Function}
    function RowJacobianMatrix(matrix::AbstractMatrix)
        S1, S2 = size(matrix)
        @assert S2 == 1 "A RowJacobianMatrix should have one column."
        new{S1}(SMatrix{S1, 1, Function}(matrix))
    end
end

RowJacobianMatrix(vector::AbstractVector) = RowJacobianMatrix(reshape(vector, :, 1))

# ColJacobianMatrix

"""
    Jacobian storing the N functions ∂f/∂x_i for f: R^N → R^M.
"""
struct ColJacobianMatrix{S} <: JacobianMatrix
    matrix::SMatrix{1, S, Function}
    function ColJacobianMatrix(matrix::AbstractMatrix)
        S1, S2 = size(matrix)
        @assert S1 == 1 "A ColJacobianMatrix should have one row."
        new{S2}(SMatrix{1, S2, Function}(matrix))
    end
end

ColJacobianMatrix(vector::AbstractVector) = ColJacobianMatrix(reshape(vector, 1, :))

# Get the partial derivative ∂f_i/∂x_j

Base.getindex(MJ::MatJacobianMatrix, i::Int, j::Int) = MJ.matrix[i, j]
Base.getindex(RJ::RowJacobianMatrix, i::Int, j::Int) = (x...) -> RJ.matrix[i](x...)[j]
Base.getindex(CJ::ColJacobianMatrix, i::Int, j::Int) = (x...) -> CJ.matrix[j](x...)[i]

# Get the gradient ∇f_i

get_gradient(MJ::MatJacobianMatrix{S1, S2}, i::Int) where {S1, S2} = (x...) -> [MJ.matrix[i, j](x...) for j in 1:S2]
get_gradient(RJ::RowJacobianMatrix, i::Int) = RJ.matrix[i]
get_gradient(CJ::ColJacobianMatrix{S}, i::Int) where {S} = (x...) -> [CJ.matrix[j](x...)[i] for j in 1:S]

# Inflate the Jacobian

inflate(MJ::MatJacobianMatrix, indices) = MatJacobianMatrix(stack(inflate(eachrow(MJ.matrix), indices), dims=1))
inflate(RJ::RowJacobianMatrix, indices) = RowJacobianMatrix(stack(inflate(eachrow(RJ.matrix), indices), dims=1))
inflate(CJ::ColJacobianMatrix, indices) = ColJacobianMatrix(map(f -> ((x...) -> inflate(f(x...), indices)), CJ.matrix))

# Append constant matrix

function append_constant(MJ::MatJacobianMatrix{S1, S2}, Z::AbstractMatrix{T}) where {S1, S2, T <: Real}
    n, m = size(Z)
    @assert S1 == n == m "The first dimension of the Jacobian matrix $(S1) must match the size of the constant matrix ($(n), $(m))."
    matrix = Matrix{Function}(undef, S1, S2 + S1)
    for i in 1:S1
        for j in 1:S2
            matrix[i, j] = (x...) -> MJ.matrix[i, j](x[1:S2]...)
        end
        for j in 1:S1
            matrix[i, S2 + j] = (x...) -> Z[i, j]
        end
    end
    return MatJacobianMatrix(matrix)
end
function append_constant(RJ::RowJacobianMatrix{S}, Z::AbstractMatrix{T}) where {S, T <: Real}
    n, m = size(Z)
    @assert S == n == m "The size of the Jacobian matrix $(S) must match the size of the constant matrix ($(n), $(m))."
    matrix = Matrix{Function}(undef, S, 1)
    for i in 1:S
        matrix[i, 1] = (x...) -> [RJ.matrix[i, 1](x...); Z[i, :]]
    end
    return RowJacobianMatrix(matrix)
end
function append_constant(CJ::ColJacobianMatrix{S}, Z::AbstractMatrix{T}) where {S, T <: Real}
    _, m = size(Z)
    # @assert
    matrix = Matrix{Function}(undef, 1, S + m)
    for j in 1:S
        matrix[1, j] = (x...) -> CJ.matrix[1, j](x[1:S]...)
    end
    for j in 1:m
        matrix[1, S + j] = (x...) -> Z[:, j]
    end
    return ColJacobianMatrix(matrix)
end