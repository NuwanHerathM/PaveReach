using LazySets, IntervalArithmetic, LinearAlgebra
using NeuralVerification

"""
    reluder(x)

    Compute the derivative of the ReLu function.
"""
function reluder(x)
    if x > 0
        return 1
    elseif x < 0
        return 0
    else
        throw(DomainError("The input to reluder must be non-zero."))
    end
end

"""
    reluder(x::IntervalArithmetic.Interval)

    For [0, 0], return [1, 1].
"""
function reluder(x::IntervalArithmetic.Interval{T}) where T <: Real
    if x.hi > 0
        upper = 1
    else 
        upper = 0
    end
    if x.lo >= 0
        lower = 1
    else 
        lower = 0
    end
    return lower <= upper ? interval(lower, upper) : interval(1, 1)
end

function act_gradient(act::NeuralVerification.ReLU, vector::AbstractVector)
    return reluder.(vector)
end

function sigmoid(x)
   return 1.0/(1.0+exp(-x))
end

function sigmoidder(x)
   return sigmoid(x)*(1.0-sigmoid(x))
end

# Needs correct rounding
function sigmoidder(x::IntervalArithmetic.Interval{T}) where T <: Real
    if x.hi < 0
        return interval(prevfloat(sigmoidder(x.lo)), nextfloat(sigmoidder(x.hi)))
    elseif x.lo > 0
        return interval(prevfloat(sigmoidder(x.hi)), nextfloat(sigmoidder(x.lo)))
    else
        low = prevfloat(min(sigmoidder(x.lo), sigmoidder(x.hi)))
        high = 0.25
        return interval(low, high)
    end
end

function act_gradient(act::NeuralVerification.Sigmoid, vector::AbstractVector)
   return sigmoidder.(vector)
end

function act_gradient(act::NeuralVerification.Id, vector::AbstractVector{T}) where T <: Real
    return ones(T, length(vector))
end

function affine_map(layer::NeuralVerification.Layer, z::AbstractVector)
    return layer.weights * z + layer.bias
end

function get_gradient(nnet::Network, x::AbstractVector)
    z = x
    gradient = Matrix(1.0LinearAlgebra.I, length(x), length(x))
    for layer in nnet.layers
        z_hat = affine_map(layer, z)
        m_gradient = act_gradient(layer.activation, z_hat)
        gradient = Diagonal(m_gradient) * layer.weights * gradient
        z = layer.activation(z_hat)
    end
    return gradient
end

function get_gradient(nnet::Network, box::IntervalBox)
    return get_gradient(nnet, box.v)
end

function has_flat_content(a::T) where T <: AbstractArray
    dimensions = size(a)
    @assert all(x -> x > 0, dimensions) "Array must have positive dimensions."
    return count(x -> x != 1, dimensions) ∈ [0, 1]
end