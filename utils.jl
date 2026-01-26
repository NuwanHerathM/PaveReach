using LazySets, IntervalArithmetic, LinearAlgebra
using NeuralVerification
# using Interpolations

"""
    act_gradient(act::NeuralVerification.ReLU, box::IntervalBox)

Compute the gradient of the ReLU activation over the given interval box.

The gradient is 1 when x > 0 and 0 when x <= 0. See the examples for details.

# Examples
```jldoctest
julia> using IntervalArithmetic, NeuralVerification
julia> box_1 = IntervalBox(interval(-2, -1))
[-2.0, -1.0]¹
julia> act_gradient(NeuralVerification.ReLU(), box_1)
1-element IntervalBox{Interval{Float64},1}:
 [0, 0]¹
julia> box_2 = IntervalBox(interval(1, 2))
[1.0, 2.0]¹
julia> act_gradient(NeuralVerification.ReLU(), box_2)
1-element IntervalBox{Interval{Float64},1}:
 [1, 1]¹
julia> box_3 = IntervalBox(interval(-1, 1))
[-1.0, 1.0]¹
julia> act_gradient(NeuralVerification.ReLU(), box_3)
1-element IntervalBox{Interval{Float64},1}:
 [0, 1]¹
julia> box_4 = IntervalBox(interval(-1, 0))
[-1.0, 0.0]¹
julia> act_gradient(NeuralVerification.ReLU(), box_4)
1-element IntervalBox{Interval{Float64},1}:
 [0, 0]¹
julia> box_5 = IntervalBox(interval(0, 1))
[0.0, 1.0]¹
julia> act_gradient(NeuralVerification.ReLU(), box_5)
1-element IntervalBox{Interval{Float64},1}:
 [0, 1]¹
```
"""
function act_gradient(act::NeuralVerification.ReLU, box::IntervalBox)
    intervals = []
    for i in 1:length(box)
        current_interval = box[i]
        if current_interval.hi > 0
            upper = 1
        else 
            upper = 0
        end
        if current_interval.lo > 0
            lower = 1
        else 
            lower = 0
        end
        push!(intervals, interval(lower, upper))
    end
    return IntervalBox(intervals)
end

function sigmoid(x)
   return 1.0/(1.0+exp(-x))
end

function sigmoid(box::IntervalBox)
   l = []
   for x in box
      push!(l, sigmoid(x))
   end
   return IntervalBox(l)
end

function sigmoidder(x::Float64)
   return sigmoid(x)*(1.0-sigmoid(x))
end

# Needs correct rounding
function sigmoidder(x::IntervalArithmetic.Interval{Float64})
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

function act_gradient(act::NeuralVerification.Sigmoid, box::IntervalBox)
   l = []
   for x in box
      push!(l, sigmoidder(x))
   end
   return IntervalBox(l)
end

act_gradient(act::NeuralVerification.Id, box::IntervalBox) = IntervalBox(interval(1), length(box))

function Diagonal end

function Diagonal(box::IntervalBox)
    M = Matrix(1.0LinearAlgebra.I, length(box), length(box))*interval(1.0,1.0)
    for i in 1:length(box)
        M[i,i] = box[i]
    end
    return M
end

function affine_map(layer::NeuralVerification.Layer, z::IntervalBox)
    return layer.weights * z + IntervalBox(layer.bias)
end

function get_gradient(nnet::Network, x::IntervalBox)
    z = x
    gradient = Matrix(1.0LinearAlgebra.I, length(x), length(x))
    for (i, layer) in enumerate(nnet.layers)
        z_hat = affine_map(layer, z)
        m_gradient = act_gradient(layer.activation, z_hat)
        gradient = Diagonal(m_gradient) * layer.weights * gradient
        z = layer.activation(z_hat)
    end
    return gradient
end

function get_gradient(nnet::Network, x::Vector{IntervalArithmetic.Interval{Float64}})
    return get_gradient(nnet, IntervalBox(x))
end