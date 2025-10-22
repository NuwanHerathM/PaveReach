#!/usr/local/bin/julia

using BenchmarkTools
using IntervalArithmetic
using Match

using Plots; pythonplot()

include("cell.jl")
include("quantifierproblem.jl")

"""
    create_is_in(qe, intervals)
# Arguments
- `qe` quantier elimination problem QuantifierProblem
- `intervals` = P_1, P_2, ..., Z
"""
function create_is_in_1(qe::QuantifierProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Forall, i) for i in 1:length(X)]..., qe.qvs...]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        R_inner, _ = QEapprox_o0(qe.f, qe.Df, dirty_quantifiers, [dirty_quantifiers for i=1:qe.n], qe.p, qe.n, [X.v..., intervals...])
        return interval(0, 0) ⊆ R_inner[1]
    end
end

function create_is_in_2(qe::QuantifierProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}}, pseudo_infinity::Int=1000, ϵ::Float64=0.1)::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Exists, i) for i in 1:length(X)]..., negation.(qe.qvs[1:end-1])..., (Exists, qe.p)]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        Z_minus = interval(-pseudo_infinity, intervals[end].lo - ϵ)
        Z_plus = interval(intervals[end].hi + ϵ, pseudo_infinity)
        _, R_outer_minus = QEapprox_o0(qe.f, qe.Df, dirty_quantifiers, [dirty_quantifiers for i=1:qe.n], qe.p, qe.n, [X.v..., intervals[1:end-1]..., Z_minus])
        _, R_outer_plus = QEapprox_o0(qe.f, qe.Df, dirty_quantifiers, [dirty_quantifiers for i=1:qe.n], qe.p, qe.n, [X.v..., intervals[1:end-1]..., Z_plus])
        return interval(0, 0) ⊈ R_outer_minus[1] && interval(0, 0) ⊈ R_outer_plus[1]
    end
end

function create_is_in_1(qe::QuantifierProblem, box::IntervalBox)::Function
    return create_is_in_1(qe, box.v)
end

function create_is_in_2(qe::QuantifierProblem, box::IntervalBox)::Function
    return create_is_in_2(qe, box.v)
end

"""
    create_is_out(qe, intervals, [pseudo_infinity, ϵ])
# Arguments
- `qe` quantifier elimination problem QuantifierProblem
- `intervals` = P_1, P_2, ..., Z
- `pseudo_infinity` a large number to represent infinity in the intervals
- `ϵ` a small number to create a margin around the interval Z for its complement
"""
function create_is_out_1(qe::QuantifierProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Exists, i) for i in 1:length(X)]..., qe.qvs...]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        _, R_outer = QEapprox_o0(qe.f, qe.Df, dirty_quantifiers, [dirty_quantifiers for i=1:qe.n], qe.p, qe.n, [X.v..., intervals...])
        return  interval(0,0) ⊈ R_outer[1]
    end
end

function create_is_out_2(qe::QuantifierProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}}, pseudo_infinity::Int=1000, ϵ::Float64=0.1)::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Forall, i) for i in 1:length(X)]..., negation.(qe.qvs[1:end-1])..., (Exists, qe.p)]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        Z_minus = interval(-pseudo_infinity, intervals[end].lo - ϵ)
        Z_plus = interval(intervals[end].hi + ϵ, pseudo_infinity)
        R_inner_minus, _ = QEapprox_o0(qe.f, qe.Df, dirty_quantifiers, [dirty_quantifiers for i=1:qe.n], qe.p, qe.n, [X.v..., intervals[1:end-1]..., Z_minus])
        R_inner_plus, _ = QEapprox_o0(qe.f, qe.Df, dirty_quantifiers, [dirty_quantifiers for i=1:qe.n], qe.p, qe.n, [X.v..., intervals[1:end-1]..., Z_plus])
        return interval(0, 0) ⊆ R_inner_minus[1] || interval(0, 0) ⊆ R_inner_plus[1]
    end
end

function create_is_out_1(qe::QuantifierProblem, box::IntervalBox)::Function
    return create_is_out_1(qe, box.v)
end

function create_is_out_2(qe::QuantifierProblem, box::IntervalBox)::Function
    return create_is_out_2(qe, box.v)
end

function is_member(cell::AbstractCell, X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
    @match typeof(cell) begin
        $CellEnd      => return cell.is_member(X)
        $AndCellStart => return all(child -> is_member(child, X), children(cell))
        $AndCell      => return all(child -> is_member(child, X), children(cell))
        _             => return any(child -> is_member(child, X), children(cell))
    end
end

function get_in_constructors_1(quantifiers)
    constructors_in = []
    is_first_forall = quantifiers[begin] == Forall
    if is_first_forall
        push!(constructors_in, AndCellStart)
        quantifiers = quantifiers[2:end]
        while quantifiers[begin] == Forall
            push!(constructors_in, AndCell)
            quantifiers = quantifiers[2:end]
        end
    else
        push!(constructors_in, OrCellStart)
        quantifiers = quantifiers[2:end]
    end
    for quantifier in quantifiers
        if quantifier == Exists
            push!(constructors_in, OrCell)
        else
            push!(constructors_in, AtomicCell)
        end
    end
    push!(constructors_in, CellEnd)
    return constructors_in
end

function get_out_constructors_1(quantifiers)
    constructors_out = []
    is_first_forall = quantifiers[begin] == Forall
    if is_first_forall
        push!(constructors_out, OrCellStart)
        quantifiers = quantifiers[2:end]
        while quantifiers[begin] == Forall
            push!(constructors_out, OrCell)
            quantifiers = quantifiers[2:end]
        end
    else
        push!(constructors_out, AndCellStart)
        quantifiers = quantifiers[2:end]
    end
    for quantifier in quantifiers
        if quantifier == Exists
            push!(constructors_out, AndCell)
        else
            push!(constructors_out, AtomicCell)
        end
    end
    push!(constructors_out, CellEnd)
    return constructors_out
end

get_in_constructors_2 = get_out_constructors_1
get_out_constructors_2 = get_in_constructors_1

function make_pz_tree_from_constructors(constructors, intervals, is_member)
    root = constructors[begin]()
    child = make_pz_tree_from_constructors(intervals, root, constructors[2:end], is_member)
    add_child(root, child)
    return root
end

function make_pz_tree_from_constructors(intervals, parent::AbstractCell, constructors, is_member)
    if length(intervals) == 1
        return constructors[begin](intervals[begin], parent, is_member)
    else
        cell = constructors[begin](intervals[begin], parent)
        child = make_pz_tree_from_constructors(intervals[2:end], cell, constructors[2:end], is_member)
        add_child(cell, child)
        return cell
    end
end

function make_in_pz_1(intervals, qe)
    is_in = create_is_in_1(qe, intervals)

    pos = [i for (_, i) in qe.qvs]
    real_pos = pos .- (qe.p - length(intervals))
    permuted_intervals = intervals[real_pos]

    quantifiers = quantifier.(qe.qvs)

    constructors_in = get_in_constructors_1(quantifiers)

    root_in = make_pz_tree_from_constructors(constructors_in, permuted_intervals, is_in)
    return root_in
end

function make_in_pz_2(intervals, qe)
    is_in = create_is_in_2(qe, intervals)

    pos = [i for (_, i) in qe.qvs]
    real_pos = pos .- (qe.p - length(intervals))
    permuted_intervals = intervals[real_pos]

    quantifiers = [negation.(quantifier.(qe.qvs[1:end-1]))..., Exists]

    constructors_in = get_in_constructors_2(quantifiers)

    root_in = make_pz_tree_from_constructors(constructors_in, permuted_intervals, is_in)
    return root_in
end

function make_out_pz_1(intervals, qe)
    is_out = create_is_out_1(qe, intervals)

    pos = [i for (_, i) in qe.qvs]
    real_pos = pos .- (qe.p - length(intervals))
    permuted_intervals = intervals[real_pos]

    quantifiers = quantifier.(qe.qvs)

    constructors_out = get_out_constructors_1(quantifiers)

    root_out = make_pz_tree_from_constructors(constructors_out, permuted_intervals, is_out)
    return root_out
end

function make_out_pz_2(intervals, qe)
    is_out = create_is_out_2(qe, intervals)

    pos = [i for (_, i) in qe.qvs]
    real_pos = pos .- (qe.p - length(intervals))
    permuted_intervals = intervals[real_pos]

    quantifiers = [negation.(quantifier.(qe.qvs[1:end-1]))..., Exists]

    constructors_out = get_out_constructors_2(quantifiers)

    root_out = make_pz_tree_from_constructors(constructors_out, permuted_intervals, is_out)
    return root_out
end

function make_pz_11(intervals, qe)
    p_in = make_in_pz_1(intervals, qe)
    p_out = make_out_pz_1(intervals, qe)
    return (p_in, p_out)
end

function make_pz_12(intervals, qe)
    p_in = make_in_pz_1(intervals, qe)
    p_out = make_out_pz_2(intervals, qe)
    return (p_in, p_out)
end

function make_pz_21(intervals, qe)
    p_in = make_in_pz_2(intervals, qe)
    p_out = make_out_pz_1(intervals, qe)
    return (p_in, p_out)
end

function make_pz_22(intervals, qe)
    p_in = make_in_pz_2(intervals, qe)
    p_out = make_out_pz_2(intervals, qe)
    return (p_in, p_out)
end

# function make_in_paving(intervals, qe)
#     is_in = create_is_in(qe, intervals)

#     pos = [i for (_, i) in qe.qvs]
#     real_pos = pos .- 1
#     permuted_intervals = intervals[real_pos]

#     quantifiers = quantifier.(qe.qvs)
#     return make_paving(permuted_intervals, quantifiers, is_in)
# end

# function make_out_paving(intervals, qe)
#     is_out = create_is_out(qe, intervals)

#     pos = [i for (_, i) in qe.qvs]
#     real_pos = pos .- 1
#     permuted_intervals = intervals[real_pos]

#     quantifiers = [negation.(quantifier.(qe.qvs[1:end-1]))..., Exists]
#     return make_paving(permuted_intervals, quantifiers, is_out)
# end

# function make_out_paving(intervals, qe)
#     is_out = create_is_out(qe, intervals)

#     pos = [i for (_, i) in qe.qvs]
#     real_pos = pos .- 1
#     permuted_intervals = intervals[real_pos]

#     quantifiers = quantifier.(qe.qvs)
#     return make_paving(permuted_intervals, quantifiers, is_out)
# end

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
function pave(p_in_0::CellStart, p_out_0::CellStart, qe::QuantifierProblem, X_0::IntervalArithmetic.IntervalBox{N, T}, ϵ::Float64, bisect_in!, bisect_out!, ratio::Float64=0.2, precision_factor::Int=10)::Tuple{Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}} where {N, T<:Number}
    inn = []
    out = []
    delta = []
    list = [(p_in_0, p_out_0, X_0)]
    while !isempty(list)
        p_in, p_out, X = pop!(list)
        if is_member(p_in, X)
            push!(inn, X)
        elseif is_member(p_out, X)
            push!(out, X)
        elseif IntervalArithmetic.diam(X) < ϵ
            push!(delta, X)
        else
            max_diam = maximum(diam, [p_in, p_out])
            if IntervalArithmetic.diam(X) > ratio * max_diam || max_diam < precision_factor * ϵ
                X_1, X_2 = bisect(X, 0.5)
                p_in_1 = deepcopy(p_in_0)
                p_out_1 = deepcopy(p_out_0)
                p_in_2 = deepcopy(p_in_0)
                p_out_2 = deepcopy(p_out_0)
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

function pave_11(p_in_0::CellStart, p_out_0::CellStart, qe::QuantifierProblem, X_0::IntervalArithmetic.IntervalBox{N, T}, ϵ::Float64, ratio::Float64=0.2, precision_factor::Int=10)::Tuple{Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}} where {N, T<:Number}
    return pave(p_in_0, p_out_0, qe, X_0, ϵ, bisect_in_1!, bisect_out_1!, ratio, precision_factor)
end

function pave_12(p_in_0::CellStart, p_out_0::CellStart, qe::QuantifierProblem, X_0::IntervalArithmetic.IntervalBox{N, T}, ϵ::Float64, ratio::Float64=0.2, precision_factor::Int=10)::Tuple{Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}} where {N, T<:Number}
    return pave(p_in_0, p_out_0, qe, X_0, ϵ, bisect_in_1!, bisect_out_2!, ratio, precision_factor)
end

function pave_21(p_in_0::CellStart, p_out_0::CellStart, qe::QuantifierProblem, X_0::IntervalArithmetic.IntervalBox{N, T}, ϵ::Float64, ratio::Float64=0.2, precision_factor::Int=10)::Tuple{Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}} where {N, T<:Number}
    return pave(p_in_0, p_out_0, qe, X_0, ϵ, bisect_in_2!, bisect_out_1!, ratio, precision_factor)
end

function pave_22(p_in_0::CellStart, p_out_0::CellStart, qe::QuantifierProblem, X_0::IntervalArithmetic.IntervalBox{N, T}, ϵ::Float64, ratio::Float64=0.2, precision_factor::Int=10)::Tuple{Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}} where {N, T<:Number}
    return pave(p_in_0, p_out_0, qe, X_0, ϵ, bisect_in_2!, bisect_out_2!, ratio, precision_factor)
end

function bisect!(root::CellStart, qe, create_is_member)
    candidate = largest_bisectable_cell(root)

    idx = depth(candidate)
    interval_1, interval_2 = bisect(candidate.interval, 0.5)

    initial_intervals = starting_intervals(root)
    intervals_1 = [i == idx ? interval_1 : initial_intervals[i] for i in 1:length(qe.qvs)]
    intervals_2 = [i == idx ? interval_2 : initial_intervals[i] for i in 1:length(qe.qvs)]

    pos = index.(qe.qvs)
    pseudo_pos = pos .- (qe.p - length(intervals))
    reordered_intervals_1 = deepcopy(initial_intervals)[pseudo_pos]
    reordered_intervals_2 = deepcopy(initial_intervals)[pseudo_pos]
    is_member_1 = create_is_member(qe, IntervalBox(reordered_intervals_1))
    is_member_2 = create_is_member(qe, IntervalBox(reordered_intervals_2))

    push!(candidate.parent, intervals_1, idx, get_constructors(root), is_member_1)
    push!(candidate.parent, intervals_2, idx, get_constructors(root), is_member_2)
    remove!(candidate)
end

function bisect_in_1!(cell::CellStart, qe::QuantifierProblem)
    return bisect!(cell, qe, create_is_in_1)
end

function bisect_in_2!(cell::CellStart, qe::QuantifierProblem)
    return bisect!(cell, qe, create_is_in_2)
end

function bisect_out_1!(cell::CellStart, qe::QuantifierProblem)
    return bisect!(cell, qe, create_is_out_1)
end

function bisect_out_2!(cell::CellStart, qe::QuantifierProblem)
    return bisect!(cell, qe, create_is_out_2)
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

function print_inn_out_delta(inn::Vector{IntervalArithmetic.IntervalBox{N, T}}, out::Vector{IntervalArithmetic.IntervalBox{N, T}}, delta::Vector{IntervalArithmetic.IntervalBox{N, T}}) where {N, T<:Number}
    if N == 1
        Base.println("Union of elements of inn: ", merge_intervals([box[1] for box in inn]))
        Base.println("Union of elements of out: ", merge_intervals([box[1] for box in out]))
        Base.println("Union of elements of delta: ", merge_intervals([box[1] for box in delta]))
    else
        Base.println("Number of elements of inn: ", length(inn))
        Base.println("Number of elements of out: ", length(out))
        Base.println("Number of elements of delta: ", length(delta))
    end
end

function draw_lines(intervals, color)
    ys = [0, 0]
    for interval in intervals
        xs = [interval.lo, interval.hi]
        plot!(xs, ys, linewidth=2, color=color, legend=:false)
    end
end


function draw_rectangles(boxes, color)
    rectangle(p, q) = Shape([p[1],q[1],q[1],p[1]], [p[2],p[2],q[2],q[2]])
    for box in boxes
        x_interval = box[1]
        y_interval = box[2]
        p = (x_interval.lo, y_interval.lo)
        q = (x_interval.hi, y_interval.hi)
        plot!(rectangle(p, q), color=color, legend=:false)
    end
end

function display(X_0, inn, out, delta)
    if isa(X_0, IntervalBox{1, <:Number})
        ylims!((-0.1,0.1))
        xticks!((X_0[1].lo:1:X_0[1].hi))

        draw_delta_lines(delta) = draw_lines(merge_intervals([box[1] for box in delta]), :yellow)
        draw_inn_lines(inn) = draw_lines(merge_intervals([box[1] for box in inn]), :green)
        draw_out_lines(out) = draw_lines(merge_intervals([box[1] for box in out]), :cyan)

        draw_delta_lines(delta)
        draw_inn_lines(inn)
        draw_out_lines(out)
        gui()
    elseif isa(X_0, IntervalBox{2, <:Number})
        xs = X_0[1]
        ys = X_0[2]
        xlims!((xs.lo, xs.hi))
        ylims!((ys.lo, ys.hi))

        draw_delta_rectangles(delta) = draw_rectangles(delta, :yellow)
        draw_inn_rectangles(inn) = draw_rectangles(inn, :green)
        draw_out_rectangles(out) = draw_rectangles(out, :cyan)

        draw_delta_rectangles(delta)
        draw_inn_rectangles(inn)
        draw_out_rectangles(out)
        gui()
    else
        error("Plotting is only supported for 1D and 2D problems.")
    end
end