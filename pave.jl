#!/usr/local/bin/julia

using BenchmarkTools
using IntervalArithmetic
using Match

include("cell.jl")
include("quantifierproblem.jl")

"""
    create_is_in(qe, intervals)
# Arguments
- `qe` quantier elimination problem QuantifierProblem
- `intervals` = P_1, P_2, ..., Z
"""
function create_is_in(qe::QuantifierProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.Interval{T}) where {T<:Number}
        quantifiers = [(Forall, 1), qe.quantifiers[1:end-1]..., (Exists, qe.p)]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        R_inner, _ = QEapprox_o0(qe.fun, dirty_quantifiers, [dirty_quantifiers for i=1:qe.n], qe.p, qe.n, [X, intervals...])
        return R_inner[1] ⊇ interval(0, 0)
    end
end

function create_is_in(qe::QuantifierProblem, box::IntervalBox)::Function
    return create_is_in(qe, box.v)
end

"""
    create_is_out(qe, intervals)
# Arguments
- `qe` quantifier elimination problem QuantifierProblem
- `intervals` = P_1, P_2, ..., Z
"""
function create_is_out(qe::QuantifierProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    pseudo_infinity = 1000
    eps = 0.1
    return function(X::IntervalArithmetic.Interval{T}) where {T<:Number}
        quantifiers = [(Forall, 1), negation.(qe.quantifiers[1:end-1])..., (Exists, qe.p)]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        Z_minus = interval(-pseudo_infinity, intervals[end].lo - eps)
        Z_plus = interval(intervals[end].hi + eps, pseudo_infinity)
        R_inner_minus, _ = QEapprox_o0(qe.fun, dirty_quantifiers, [dirty_quantifiers for i=1:qe.n], qe.p, qe.n, [X, intervals[1:end-1]..., Z_minus])
        R_inner_plus, _ = QEapprox_o0(qe.fun, dirty_quantifiers, [dirty_quantifiers for i=1:qe.n], qe.p, qe.n, [X, intervals[1:end-1]..., Z_plus])
        return  R_inner_minus[1] ⊇ interval(0, 0) || R_inner_plus[1] ⊇ interval(0, 0)
    end
end

function create_is_out(qe::QuantifierProblem, box::IntervalBox)::Function
    return create_is_out(qe, box.v)
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

function pave(p_0::MembershipCell, qe::QuantifierProblem, X_0::IntervalArithmetic.Interval{T}, ϵ::Float64)::Tuple{Vector{IntervalArithmetic.Interval{T}}, Vector{IntervalArithmetic.Interval{T}}, Vector{IntervalArithmetic.Interval{T}}} where {T<:Number}
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
        elseif IntervalArithmetic.diam(X) < ϵ
            push!(delta, X)
        else
            if IntervalArithmetic.diam(X) > diam(p) || diam(p) < 10 * ϵ
                X_1, X_2 = bisect(X, 0.5)
                push!(list, (p, X_1))
                push!(list, (p, X_2))
            else
                bisect!(p, qe)
                push!(list, (p, X))
            end
        end
    end
    return (inn, out, delta)
end

function bisect!(cell::MembershipCell, qe::QuantifierProblem)
    if isleaf(cell)
        box_1, box_2 = bisect(cell.box, 0.5)
        
        is_in_1 = create_is_in(qe, box_1)
        is_out_1 = create_is_out(qe, box_1)
        cell_1 = make_membershipcell_leaf(box_1, is_in_1, is_out_1, cell)
        
        is_in_2 = create_is_in(qe, box_2)
        is_out_2 = create_is_out(qe, box_2)
        cell_2 = make_membershipcell_leaf(box_2, is_in_2, is_out_2, cell)
        
        dim = 1
        while box_1[dim] == box_2[dim]
            dim += 1
        end
        dim_X = 1
        q = quantifier(qe, dim + dim_X)
        if q == Exists
            cell.is_in = X -> cell_1.is_in(X) || cell_2.is_in(X)
            cell.is_out = X -> cell_1.is_out(X) && cell_2.is_out(X)
        elseif q == Forall
            cell.is_in = X -> cell_1.is_in(X) && cell_2.is_in(X)
            cell.is_out = X -> cell_1.is_out(X) || cell_2.is_out(X)
        else
            error("Unknown quantifier: $q")
        end
    else
        heights = height.(cell.children)
        if heights[1] > heights[2]
            bisect!(cell.children[2], qe)
        else
            bisect!(cell.children[1], qe)
        end
    end
end

function pave(qe::QuantifierProblem, intervals::Vector{IntervalArithmetic.Interval{T}}, X_0::IntervalArithmetic.Interval{T}, ϵ::Float64) where {T<:Number}
    is_in = create_is_in(qe, intervals)
    is_out = create_is_out(qe, intervals)
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
    Base.println("Union of elements of inn: ", merge_intervals(inn))
    Base.println("Union of elements of out: ", merge_intervals(out))
    Base.println("Union of elements of delta: ", merge_intervals(delta))
end