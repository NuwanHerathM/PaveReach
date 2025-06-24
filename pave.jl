#!/usr/local/bin/julia

using BenchmarkTools
using IntervalArithmetic

include("cell.jl")

function negation(quantifiers)
    negated_quantifiers = []
    for i in eachindex(quantifiers)
        if quantifiers[i] == "forall"
            push!(negated_quantifiers, "exists")
        elseif quantifiers[i] == "exists"
            push!(negated_quantifiers, "forall")
        else
            push!(negated_quantifiers, quantifiers[i])
        end
    end
    return negated_quantifiers
end

"""
    create_is_in(QE, intervals)
# Arguments
- `QE` = (g, quantifiers, p, n) where:
  - `g` is the formula to be evaluated
  - `quantifiers` is a list of quantifiers in the formula
  - `p` is the number of variables in the formula
  - `n` is the number of intervals
- `intervals` = P_1, P_2, ..., Z
"""
function create_is_in(QE, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.Interval{T}) where {T<:Number}
        quantifiers = ["forall", 1, QE[2][1:end-2]..., "exists", QE[3]]
        R_inner, _ = QEapprox_o0(QE[1], quantifiers, [quantifiers for i=1:QE[4]], QE[3], QE[4], [X, intervals...])
        return R_inner[1] ⊇ interval(0, 0)
    end
end

function create_is_in(QE, box::IntervalBox)::Function
    return create_is_in(QE, box.v)
end

"""
    create_is_out(QE, intervals)
# Arguments
- `QE` = (g, quantifiers, p, n) where:
  - `g` is the formula to be evaluated
  - `quantifiers` is a list of quantifiers in the formula
  - `p` is the number of variables in the formula
  - `n` is the number of intervals
- `intervals` = P_1, P_2, ..., Z
"""
function create_is_out(QE, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    pseudo_infinity = 1000
    eps = 0.1
    return function(X::IntervalArithmetic.Interval{T}) where {T<:Number}
        quantifiers = ["forall", 1, negation(QE[2][1:end-2])..., "exists", QE[3]]
        Z_minus = interval(-pseudo_infinity, intervals[end].lo - eps)
        Z_plus = interval(intervals[end].hi + eps, pseudo_infinity)
        R_inner_minus, _ = QEapprox_o0(QE[1], quantifiers, [quantifiers for i=1:QE[4]], QE[3], QE[4], [X, intervals[1:end-1]..., Z_minus])
        R_inner_plus, _ = QEapprox_o0(QE[1], quantifiers, [quantifiers for i=1:QE[4]], QE[3], QE[4], [X, intervals[1:end-1]..., Z_plus])
        return  R_inner_minus[1] ⊇ interval(0, 0) || R_inner_plus[1] ⊇ interval(0, 0)
    end
end

function create_is_out(QE, box::IntervalBox)::Function
    return create_is_out(QE, box.v)
end

function pave(is_in::Function, is_out::Function, X_0::IntervalArithmetic.Interval{T}, ϵ::Float64)::Tuple{Vector{IntervalArithmetic.Interval{T}}, Vector{IntervalArithmetic.Interval{T}}, Vector{IntervalArithmetic.Interval{T}}} where {T<:Number}
    inn = []
    out = []
    delta = []
    list = [X_0]
    while !isempty(list)
        X = pop!(list)
        if is_in(X)
            push!(inn, X)
        elseif is_out(X)
            push!(out, X)
        elseif diam(X) < ϵ
            push!(delta, X)
        else
            X_1, X_2 = bisect(X, 0.5)
            push!(list, X_1)
            push!(list, X_2)
        end
    end
    return (inn, out, delta)
end

function pave(p_0::MembershipCell, QE, X_0::IntervalArithmetic.Interval{T}, ϵ::Float64)::Tuple{Vector{IntervalArithmetic.Interval{T}}, Vector{IntervalArithmetic.Interval{T}}, Vector{IntervalArithmetic.Interval{T}}} where {T<:Number}
    inn = []
    out = []
    delta = []
    list = [(p_0, X_0)]
    while !isempty(list)
        p, X = pop!(list)
        if p.is_in(X)
            push!(inn, X)
        elseif p.is_out(X)
            push!(out, X)
        elseif diam(X) < ϵ
            push!(delta, X)
        else
            if diam(X) > diamm(p) || diamm(p) < 10 * ϵ
                X_1, X_2 = bisect(X, 0.5)
                push!(list, (p, X_1))
                push!(list, (p, X_2))
            else
                bisect!(p, QE)
                push!(list, (p, X))
            end
        end
    end
    return (inn, out, delta)
end

# function get_quantifier_from_QE(QE::Tuple{Vector{Num}, Vector{Any}, Int64, Int64}, dim::Int)
function get_quantifier_from_QE(QE, dim::Int)
    i = 2
    while QE[2][i] != dim
        i += 2
    end
    return QE[2][i-1]
end

function bisect!(cell::MembershipCell, QE)
    if isleaf(cell)
        box_1, box_2 = bisect(cell.box, 0.5)
        
        is_in_1 = create_is_in(QE, box_1)
        is_out_1 = create_is_out(QE, box_1)
        cell_1 = make_membershipcell_leaf(box_1, is_in_1, is_out_1, cell)
        
        is_in_2 = create_is_in(QE, box_2)
        is_out_2 = create_is_out(QE, box_2)
        cell_2 = make_membershipcell_leaf(box_2, is_in_2, is_out_2, cell)
        
        dim = 1
        while box_1[dim] == box_2[dim]
            dim += 1
        end
        quantifier = get_quantifier_from_QE(QE, dim + 1)
        if quantifier == "exists"
            cell.is_in = X -> cell_1.is_in(X) || cell_2.is_in(X)
            cell.is_out = X -> cell_1.is_out(X) && cell_2.is_out(X)
        else
            cell.is_in = X -> cell_1.is_in(X) && cell_2.is_in(X)
            cell.is_out = X -> cell_1.is_out(X) || cell_2.is_out(X)
        end
    else
        heights = height.(cell.children)
        if heights[1] > heights[2]
            bisect!(cell.children[2], QE)
        else
            bisect!(cell.children[1], QE)
        end
    end
end

# function pave(QE::Tuple{Vector{Num}, Vector{Any}, Int64, Int64}, intervals::Vector{IntervalArithmetic.Interval{T}}, X_0::IntervalArithmetic.Interval{T}, ϵ::Float64) where {T<:Number}
function pave(QE, intervals::Vector{IntervalArithmetic.Interval{T}}, X_0::IntervalArithmetic.Interval{T}, ϵ::Float64) where {T<:Number}
    is_in = create_is_in(QE, intervals)
    is_out = create_is_out(QE, intervals)
    return pave(is_in, is_out, X_0, ϵ)
end

function merge_intervals(intervals::Vector{IntervalArithmetic.Interval{T}}) where {T<:Number}
    disconnected = deepcopy(intervals)
    merged = []
    while !isempty(disconnected)
        current = pop!(disconnected)
        i = 1
        while i <= length(disconnected)
            if current.lo == disconnected[i].hi
                current = interval(disconnected[i].lo, current.hi)
                popat!(disconnected, i)
                push!(disconnected, current)
                break
            elseif current.hi == disconnected[i].lo
                current = interval(current.lo, disconnected[i].hi)
                popat!(disconnected, i)
                push!(disconnected, current)
                break
            else
                i += 1
            end
        end
        if i > length(disconnected)
            push!(merged, current)
        end
    end
    return merged
end

function print_inn_out_delta(inn::Vector{IntervalArithmetic.Interval{T}}, out::Vector{IntervalArithmetic.Interval{T}}, delta::Vector{IntervalArithmetic.Interval{T}}) where {T<:Number}
    println("Union of elements of inn: ", merge_intervals(inn))
    println("Union of elements of out: ", merge_intervals(out))
    println("Union of elements of delta: ", merge_intervals(delta))
end