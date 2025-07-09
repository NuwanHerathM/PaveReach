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
        R_inner, _ = QEapprox_o0(qe.f, qe.Df, dirty_quantifiers, [dirty_quantifiers for i=1:qe.n], qe.p, qe.n, [X, intervals...])
        # println("R_inner: ", R_inner)
        return R_inner[1] ⊇ interval(0, 0)
    end
end

function create_is_in(qe::QuantifierProblem, box::IntervalBox)::Function
    return create_is_in(qe, box.v)
end

"""
    create_is_out(qe, intervals, [pseudo_infinity, ϵ])
# Arguments
- `qe` quantifier elimination problem QuantifierProblem
- `intervals` = P_1, P_2, ..., Z
- `pseudo_infinity` a large number to represent infinity in the intervals
- `ϵ` a small number to create a margin around the interval Z for its complement
"""
function create_is_out(qe::QuantifierProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}}, pseudo_infinity::Int=1000, ϵ::Float64=0.1)::Function where {T<:Number}
    return function(X::IntervalArithmetic.Interval{T}) where {T<:Number}
        quantifiers = [(Forall, 1), negation.(qe.quantifiers[1:end-1])..., (Exists, qe.p)]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        Z_minus = interval(-pseudo_infinity, intervals[end].lo - ϵ)
        Z_plus = interval(intervals[end].hi + ϵ, pseudo_infinity)
        R_inner_minus, _ = QEapprox_o0(qe.f, qe.Df, dirty_quantifiers, [dirty_quantifiers for i=1:qe.n], qe.p, qe.n, [X, intervals[1:end-1]..., Z_minus])
        R_inner_plus, _ = QEapprox_o0(qe.f, qe.Df, dirty_quantifiers, [dirty_quantifiers for i=1:qe.n], qe.p, qe.n, [X, intervals[1:end-1]..., Z_plus])
        return  R_inner_minus[1] ⊇ interval(0, 0) || R_inner_plus[1] ⊇ interval(0, 0)
    end
end

function create_is_out(qe::QuantifierProblem, box::IntervalBox)::Function
    return create_is_out(qe, box.v)
end

function is_in(cell::AbstractCell, X::IntervalArithmetic.Interval{T}, quantifiers::Vector{Tuple{Quantifier, Int}}, i::Int) where T<:Number
    if !isleaf(cell)
        q, _ = quantifiers[i]
        if q == Exists
            return any(child -> is_in(child, X, quantifiers, i+1), cell.children)
        else
            return all(child -> is_in(child, X, quantifiers, i+1), cell.children)
        end
    else
        return cell.is_in(X)
    end
end

function is_in(cell::CellStart, X::IntervalArithmetic.Interval{T}, quantifiers::Vector{Tuple{Quantifier, Int}}) where T<:Number
    return is_in(cell, X, quantifiers, 1)
end

function is_out(cell::AbstractCell, X::IntervalArithmetic.Interval{T}, quantifiers::Vector{Tuple{Quantifier, Int}}, i::Int) where T<:Number
    if !isleaf(cell)
        q, _ = quantifiers[i]
        if q == Exists
            return all(child -> is_out(child, X, quantifiers, i+1), cell.children)
        else
            return any(child -> is_out(child, X, quantifiers, i+1), cell.children)
        end
    else
        return cell.is_out(X)
    end
    
end

function is_out(cell::CellStart, X::IntervalArithmetic.Interval{T}, quantifiers::Vector{Tuple{Quantifier, Int}}) where T<:Number
    return is_out(cell, X, quantifiers, 1)
end

function make_paving(intervals, is_in, is_out)
    root = CellStart()
    child = make_paving(intervals, root, is_in, is_out)
    push!(root.children, child)
    return root
end

function make_paving(intervals, qe)
    is_in = create_is_in(qe, intervals)
    is_out = create_is_out(qe, intervals)

    pos = [i for (_, i) in qe.quantifiers]
    real_pos = pos .- 1
    permuted_intervals = intervals[real_pos]

    return make_paving(permuted_intervals, is_in, is_out)
end

function make_paving(intervals, parent::CellStart, is_in, is_out)
    return make_paving(intervals, 1, parent, is_in, is_out)
end

function make_paving(intervals, i, parent::Union{CellStart, Cell}, is_in, is_out)
    interval = intervals[i]
    if interval != intervals[end]
        current = Cell(interval, parent)
        child = make_paving(intervals, i+1, current, is_in, is_out)
        push!(current.children, child)
        return current
    else
        return CellEnd(interval, parent, is_in, is_out)
    end
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
        elseif IntervalArithmetic.diam(X) < ϵ
            push!(delta, X)
        else
            X_1, X_2 = bisect(X, 0.5)
            push!(list, X_1)
            push!(list, X_2)
        end
    end
    return (inn, out, delta)
end

"""
    pave(p_0, qe, X_0, ϵ, [ratio, precision_factor])
# Arguments
- `p_0` initial p
- `qe` quantifier elimination problem
- `X_0` initial interval
- `ϵ` precision threshold for X: if diam(X) < ϵ then X is not split
- `ratio` ratio to determine when to split the cell: if diam(X) / diam(p) <= ratio then split p
- `precision_factor` factor to get the precision threshold for p: if diam(p) < precision_factor * ϵ then p is not split
"""
function pave(p_0::CellStart, qe::QuantifierProblem, X_0::IntervalArithmetic.Interval{T}, ϵ::Float64, ratio::Float64=0.2, precision_factor::Int=10)::Tuple{Vector{IntervalArithmetic.Interval{T}}, Vector{IntervalArithmetic.Interval{T}}, Vector{IntervalArithmetic.Interval{T}}} where {T<:Number}
    inn = []
    out = []
    delta = []
    list = [(p_0, X_0)]
    while !isempty(list)
        p, X = pop!(list)
        if is_in(p, X, qe.quantifiers)
            push!(inn, X)
        elseif is_out(p, X, qe.quantifiers)
            # println(X)
            # println(width(p))
            # print_tree(p)
            push!(out, X)
            if X.lo == 3.75 && X.hi == 3.875
                return p
            end
        elseif IntervalArithmetic.diam(X) < ϵ
            push!(delta, X)
        else
            if IntervalArithmetic.diam(X) > ratio * diam(p) || diam(p) < precision_factor * ϵ
                X_1, X_2 = bisect(X, 0.5)
                p_1 = deepcopy(p_0)
                p_2 = deepcopy(p_0)
                push!(list, (p_1, X_1))
                push!(list, (p_2, X_2))
            else
                bisect!(p, qe)
                push!(list, (p, X))
            end
        end
    end
    return (inn, out, delta)
end

function bisect!(cell::CellStart, qe::QuantifierProblem)
    intervals = []

    current = cell
    while !isleaf(current)
        diams = diam.(current.children)
        pos = argmax(diams)
        current = current.children[pos]
        push!(intervals, current.interval)
    end
    
    box = IntervalBox(intervals)
    box_1, box_2 = bisect(box, 0.5)

    pos = [i for (q, i) in qe.quantifiers]
    pseudo_pos = pos .- 1
    intervals_1 = box_1.v[pseudo_pos]
    intervals_2 = box_2.v[pseudo_pos]
    is_in_1 = create_is_in(qe, intervals_1)
    is_out_1 = create_is_out(qe, intervals_1)
    is_in_2 = create_is_in(qe, intervals_2)
    is_out_2 = create_is_out(qe, intervals_2)

    remove!(cell, box.v)
    push!(cell, box_1.v, is_in_1, is_out_1)
    push!(cell, box_2.v, is_in_2, is_out_2)
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