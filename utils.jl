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