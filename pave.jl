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
        quantifiers = [(Forall, 1), qe.qvs[1:end-1]..., (Exists, qe.p)]
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
        quantifiers = [(Forall, 1), negation.(qe.qvs[1:end-1])..., (Exists, qe.p)]
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

# function is_in(cell::AbstractCell, X::IntervalArithmetic.Interval{T}, quantifiers::Vector{Quantifier}, i::Int) where T<:Number
#     if !isleaf(cell)
#         q = quantifiers[i]
#         if q == Exists
#             return any(child -> is_in(child, X, quantifiers, i+1), cell.children)
#         else
#             return all(child -> is_in(child, X, quantifiers, i+1), cell.children)
#         end
#     else
#         return cell.is_in(X)
#     end
# end

# function is_in(cell::CellStart, X::IntervalArithmetic.Interval{T}, quantifiers::Vector{Quantifier}) where T<:Number
#     return is_in(cell, X, quantifiers, 1)
# end

# function is_out(cell::AbstractCell, X::IntervalArithmetic.Interval{T}, quantifiers::Vector{Quantifier}, i::Int) where T<:Number
#     if !isleaf(cell)
#         q = quantifiers[i]
#         if q == Exists
#             return all(child -> is_out(child, X, quantifiers, i+1), cell.children)
#         else
#             return any(child -> is_out(child, X, quantifiers, i+1), cell.children)
#         end
#     else
#         return cell.is_out(X)
#     end
    
# end

# function is_out(cell::CellStart, X::IntervalArithmetic.Interval{T}, quantifiers::Vector{Quantifier}) where T<:Number
#     return is_out(cell, X, quantifiers, 1)
# end

function is_member(cell::AbstractCell, X::IntervalArithmetic.Interval{T}) where T<:Number
    if isleaf(cell)
        return cell.is_member(X)
    else
        return any(child -> is_member(child, X), cell.children)
    end
end

function make_paving(intervals, quantifiers, is_member)
    root = CellStart()
    child = make_paving(intervals, root, quantifiers, is_member)
    push!(root.children, child)
    return root
end

function make_in_paving(intervals, qe)
    is_in = create_is_in(qe, intervals)

    pos = [i for (_, i) in qe.qvs]
    real_pos = pos .- 1
    permuted_intervals = intervals[real_pos]

    quantifiers = quantifier.(qe.qvs)
    return make_paving(permuted_intervals, quantifiers, is_in)
end

function make_out_paving(intervals, qe)
    is_out = create_is_out(qe, intervals)

    pos = [i for (_, i) in qe.qvs]
    real_pos = pos .- 1
    permuted_intervals = intervals[real_pos]

    quantifiers = negation.(quantifier.(qe.qvs))
    return make_paving(permuted_intervals, quantifiers, is_out)
end

function make_paving(intervals, parent::CellStart, quantifiers, is_member)
    return make_paving(intervals, 1, parent, quantifiers, is_member)
end

function make_paving(intervals, i, parent::Union{CellStart, Cell}, quantifiers, is_member)
    interval = intervals[i]
    q = quantifiers[i]
    if interval != intervals[end]
        if q == Exists
            current = BisectableCell(interval, parent)
        else
            current = AtomicCell(interval, parent)
        end
        child = make_paving(intervals, i+1, current, quantifiers, is_member)
        push!(current.children, child)
        return current
    else
        if q == Exists
            return BisectableCellEnd(interval, parent, is_member)
        else
            return AtomicCellEnd(interval, parent, is_member)
        end
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
function pave(p_in_0::CellStart, p_out_0::CellStart, qe::QuantifierProblem, X_0::IntervalArithmetic.Interval{T}, ϵ::Float64, ratio::Float64=0.2, precision_factor::Int=10)::Tuple{Vector{IntervalArithmetic.Interval{T}}, Vector{IntervalArithmetic.Interval{T}}, Vector{IntervalArithmetic.Interval{T}}} where {T<:Number}
    inn = []
    out = []
    delta = []
    list = [(p_in_0, p_out_0, X_0)]
    while !isempty(list)
        p_in, p_out, X = pop!(list)
        if is_member(p_in, X)
            push!(inn, X)
        elseif is_member(p_out, X)
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
            # diam = maximum(diam(p_in), diam(p_out))
            if IntervalArithmetic.diam(X) > ratio * diam(p_in) || diam(p_in) < precision_factor * ϵ
                X_1, X_2 = bisect(X, 0.5)
                p_in_1 = deepcopy(p_in) # ?
                p_out_1 = deepcopy(p_out) # ?
                p_in_2 = deepcopy(p_in)
                p_out_2 = deepcopy(p_out)
                push!(list, (p_in_1, p_out_1, X_1))
                push!(list, (p_in_2, p_out_2, X_2))
            else
                bisect_in!(p_in, qe)
                bisect_out!(p_out, qe)
                push!(list, (p_in, p_out, X))
            end
        end
    end
    return (inn, out, delta)
end

# import Base: push!

function _bisect(cell::CellStart)
    n = height(cell)
    bisectable_intervals = []
    atomic_intervals = []
    bisectable_indices = []
    i = 1

    current = cell
    while !isleaf(current)
        diams = diam.(current.children)
        pos = argmax(diams)
        current = current.children[pos]
        if isa(current, Bisectable)
            push!(bisectable_intervals, current.interval)
            push!(bisectable_indices, i)
        else
            push!(atomic_intervals, current.interval)
        end
        i += 1
    end
    
    if isempty(bisectable_intervals)
        return (atomic_intervals, [], current)
    else
        subbox = IntervalBox(bisectable_intervals)
        subbox_1, subbox_2 = bisect(subbox, 0.5)

        intervals_1 = []
        intervals_2 = []
        i_bisectable = 0
        i_atomic = 0
        while i_bisectable + i_atomic < n
            if i_bisectable + i_atomic + 1 in bisectable_indices
                i_bisectable += 1
                push!(intervals_1, subbox_1.v[i_bisectable])
                push!(intervals_2, subbox_2.v[i_bisectable])
            else
                i_atomic += 1
                push!(intervals_1, atomic_intervals[i_atomic])
                push!(intervals_2, atomic_intervals[i_atomic])
            end
        end

        return (intervals_1, intervals_2, current)
    end
end

function bisect!(cell::CellStart, qe::QuantifierProblem, create_is_member::Function)
    intervals_1, intervals_2, cell_end = _bisect(cell)

    if isempty(intervals_2)
        # No bisectable intervals found, just return
        return
    end

    pos = index.(qe.qvs)
    pseudo_pos = pos .- 1
    reordered_intervals_1 = intervals_1[pseudo_pos]
    reordered_intervals_2 = intervals_2[pseudo_pos]
    is_member_1 = create_is_member(qe, IntervalBox(reordered_intervals_1))
    is_member_2 = create_is_member(qe, IntervalBox(reordered_intervals_2))

    # print_tree(cell)
    # println(intervals_1)
    push!(cell, intervals_1, is_member_1)
    push!(cell, intervals_2, is_member_2)
    remove!(cell_end)
end

function bisect_in!(cell::CellStart, qe::QuantifierProblem)
    return bisect!(cell, qe, create_is_in)
end

function bisect_out!(cell::CellStart, qe::QuantifierProblem)
    return bisect!(cell, qe, create_is_out)
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