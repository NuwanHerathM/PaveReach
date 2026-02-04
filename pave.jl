using Plots
using Luxor
using MathTeXEngine

include("genreach2.jl")
include("quantifiedconstraintproblem.jl")

# Paving

function bisect_eps(interval, ϵ)
    parts = [interval]
    while diam(first(parts)) > ϵ
        newparts = []
        for current in parts
            a, b = IntervalArithmetic.bisect(current)
            push!(newparts, a)
            push!(newparts, b)
        end
        parts = newparts
    end
    return parts
end

parts = bisect_eps(interval(0, 1), 0.1)

function bisect_eps_quantifier!(intervals, qvs, eps, p, n, quantifier)
    @assert sum(length.(intervals); init=0) ==  length(intervals) "Each interval should be a single interval."
    pos_quantifier = [i for (q, i) in qvs if q == quantifier] .- (p - n - length(qvs))

    for i in pos_quantifier
        intervals[i] = bisect_eps(intervals[i][1], eps)
    end
end

bisect_eps_exists!(intervals, qvs, eps, p, n) = bisect_eps_quantifier!(intervals, qvs, eps, p, n, Exists)
bisect_eps_forall!(intervals, qvs, eps, p, n) = bisect_eps_quantifier!(intervals, qvs, eps, p, n, Forall)

function pointify_quantifier!(intervals, qvs, n, quantifier)
    pos_quantifier = [i for (q, i) in qvs if q == quantifier] .- (p - n - length(qvs))

    for i in pos_quantifier
        intervals[i] = interval.(mid.(intervals[i]))
    end
end

pointify_exists!(intervals, qvs, n) = pointify_quantifier!(intervals, qvs, n, Exists)
pointify_forall!(intervals, qvs, n) = pointify_quantifier!(intervals, qvs, n, Forall)

function refine_in!(p_in, qvs, eps, p, n)
    bisect_eps_exists!(p_in, qvs, eps, p, n)
    pointify_exists!(p_in, qvs, n)
end

function refine_out!(p_out, qvs, eps, p, n)
    bisect_eps_forall!(p_out, qvs, eps, p,n)
    pointify_forall!(p_out, qvs, n)
end

# function bisect_largest!(intervals)
#     (_, pos_max) = findmax(IntervalArithmetic.diam.(first.(intervals)))
#     parts = []
#     for interval in intervals[pos_max]
#         a, b = IntervalArithmetic.bisect(interval)
#         push!(parts, a)
#         push!(parts, b)
#     end
#     intervals[pos_max] = parts
# end

function bisect_largest_quantifier!(intervals, qvs, p, n, quantifier)
    pos_quantifier = [i for (q, i) in qvs if q == quantifier] .- (p - n - length(qvs))
    diams = [if (i in pos_quantifier) IntervalArithmetic.diam(first(intervals[i])) else -1 end for i in 1:length(intervals)]
    (_, pos_max) = findmax(diams)
    parts = []
    for interval in intervals[pos_max]
        a, b = IntervalArithmetic.bisect(interval)
        push!(parts, a)
        push!(parts, b)
    end
    intervals[pos_max] = parts
end

bisect_largest_exists!(intervals, qvs, p, n) = bisect_largest_quantifier!(intervals, qvs, p, n, Exists)
bisect_largest_forall!(intervals, qvs, p, n) = bisect_largest_quantifier!(intervals, qvs, p, n, Forall)

function increment!(indices, lengths, pos, i)
    if length(pos) == 0
        return
    end
    if i == length(pos)
        indices[pos[begin:end-1]] = lengths[pos[begin:end-1]]
        indices[pos[end]] = lengths[pos[end]] + 1
        return
    end
    indices[pos[end-i]] += 1
    is_remainder = indices[pos[end-i]] > lengths[pos[end-i]]
    if is_remainder
        indices[pos[end-i]] = 1
        increment!(indices, lengths, pos, i+1)
    end
end

increment!(indices, lengths, pos) = increment!(indices, lengths, pos, 0)

function create_is_in_1(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Forall, i) for i in 1:length(X)]..., qcp.qvs..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        quantifiers_relaxed = [[[(Forall, i) for i in 1:length(X)]..., qcp.qvs_relaxed[j]..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...] for j in 1:qcp.n]
        dirty_qs = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        problem = qcp.problem
        R_inner, _ = QEapprox_o0(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qcp.p, qcp.n, [X.v..., intervals...])
        if any(isempty, R_inner)
            return false
        end
        return all(interval(0, 0) ⊆ interval(min(R_inner[i]), max(R_inner[i])) for i in 1:qcp.n)
    end
end

function create_is_in_2(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}}, f_bounds::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Exists, i) for i in 1:length(X)]..., negation.(qcp.qvs)..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        quantifiers_relaxed = [[[(Exists, i) for i in 1:length(X)]..., negation.(qcp.qvs_relaxed[j])..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...] for j in 1:qcp.n]
        dirty_qs = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        problem = qcp.problem
        G_minus = [interval(-∞, intervals[end-i].lo) ∩ f_bounds[end-i] for i in (qcp.n-1):-1:0]
        if any(isempty, G_minus)
            test_minus = true
        else
            _, R_outer_minus = QEapprox_o0(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qcp.p, qcp.n, [X.v..., intervals[1:end-qcp.n]..., G_minus...])
            test_minus = any([interval(0, 0) ⊈ interval(min(R_outer_minus[i]), max(R_outer_minus[i])) for i in 1:qcp.n])
        end
        G_plus = [interval(intervals[end-i].hi, ∞) ∩ f_bounds[end-i] for i in (qcp.n-1):-1:0]
        if any(isempty, G_plus)
            test_plus = true
        else
            _, R_outer_plus = QEapprox_o0(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qcp.p, qcp.n, [X.v..., intervals[1:end-qcp.n]..., G_plus...])
            test_plus = any([interval(0, 0) ⊈ interval(min(R_outer_plus[i]), max(R_outer_plus[i])) for i in 1:qcp.n])
        end
        return test_minus && test_plus
    end
end

function create_is_out_1(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Exists, i) for i in 1:length(X)]..., qcp.qvs..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        quantifiers_relaxed = [[[(Exists, i) for i in 1:length(X)]..., qcp.qvs_relaxed[j]..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...] for j in 1:qcp.n]
        dirty_qs = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        problem = qcp.problem
        _, R_outer = QEapprox_o0(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qcp.p, qcp.n, [X.v..., intervals...])
        if all(isempty, R_outer)
            return true
        end
        return  any([interval(0,0) ⊈ interval(min(R_outer[i]), max(R_outer[i])) for i in 1:qcp.n if !isempty(R_outer[i])])
    end
end

function create_is_out_2(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}},  f_bounds::AbstractVector{IntervalArithmetic.Interval{T}}, ϵ::Float64=0.1)::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Forall, i) for i in 1:length(X)]..., negation.(qcp.qvs)..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        quantifiers_relaxed = [[[(Forall, i) for i in 1:length(X)]..., negation.(qcp.qvs_relaxed[j])..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...] for j in 1:qcp.n]
        dirty_qs = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        problem = qcp.problem
        G_minus = [interval(-∞, intervals[end-i].lo - ϵ) ∩ f_bounds[end-i] for i in (qcp.n-1):-1:0]
        if any(isempty, G_minus)
            test_minus = false
        else
            R_inner_minus, _ = QEapprox_o0(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qcp.p, qcp.n, [X.v..., intervals[1:end-qcp.n]..., G_minus...])
            if any(isempty, R_inner_minus)
                test_minus = false
            else
                test_minus = all([interval(0, 0) ⊆ interval(min(R_inner_minus[i]), max(R_inner_minus[i])) for i in 1:qcp.n])
            end
        end
        G_plus = [interval(intervals[end-i].hi + ϵ, ∞) ∩ f_bounds[end-i] for i in (qcp.n-1):-1:0]
        if any(isempty, G_plus)
            test_plus = false
        else
            R_inner_plus, _ = QEapprox_o0(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qcp.p, qcp.n, [X.v..., intervals[1:end-qcp.n]..., G_plus...])
            if any(isempty, R_inner_plus)
                test_plus = false
            else
                test_plus = all([interval(0, 0) ⊆ interval(min(R_inner_plus[i]), max(R_inner_plus[i])) for i in 1:qcp.n])
            end
        end
        return test_minus || test_plus
    end
end

function bounds(f, X, interval)
    return [f[i]([X.v..., interval...]) for i in 1:length(f)]
end

function check_is_in(X_0, p_in, G, qcp, criterion)
    @assert criterion == 1 || criterion == 2

    indices_forall = [i for (q, i) in qcp.qvs if q == Forall] .- length(X_0)
    indices_exists = [i for (q, i) in qcp.qvs if q == Exists] .- length(X_0)

    indices = [1 for i in 1:length(p_in)]
    lengths = length.(p_in)
    is_in_union = false
    while !is_in_union && (isempty(indices_exists) || indices[indices_exists] <= lengths[indices_exists])
        is_in_intersection = true
        while is_in_intersection && (isempty(indices_forall) || indices[indices_forall] <= lengths[indices_forall])
            sub_interval = [[p_in[i][indices[i]] for i in 1:length(p_in)]..., G...]
            if criterion == 1
                is_in = create_is_in_1(qcp, sub_interval)
            end
            if criterion == 2
                f_bounds = bounds(qcp.problem.f, X_0, sub_interval)
                is_in = create_is_in_2(qcp, sub_interval, f_bounds)
            end
            is_in_intersection &= is_in(X_0)
            if isempty(indices_forall)
                break
            end
            increment!(indices, lengths, indices_forall)
        end
        for i in indices_forall
            indices[i] = 1
        end
        is_in_union |= is_in_intersection
        if isempty(indices_exists)
            break
        end
        increment!(indices, lengths, indices_exists)
    end
    return is_in_union
end

check_is_in_1(X_0, p_in, G, qcp) = check_is_in(X_0, p_in, G, qcp, 1)
check_is_in_2(X_0, p_in, G, qcp) = check_is_in(X_0, p_in, G, qcp, 2)

function check_is_out(X_0, p_out, G, qcp, criterion)
    @assert criterion == 1 || criterion == 2

    indices_forall = [i for (q, i) in qcp.qvs if q == Forall] .- length(X_0)
    indices_exists = [i for (q, i) in qcp.qvs if q == Exists] .- length(X_0)

    indices = [1 for i in 1:length(p_in)]
    lengths = length.(p_out)
    is_out_intersection = true
    while is_out_intersection && (isempty(indices_exists) || indices[indices_exists] <= lengths[indices_exists])
        is_out_union = false
        while !is_out_union && (isempty(indices_forall) || indices[indices_forall] <= lengths[indices_forall])
            sub_interval = [[p_out[i][indices[i]] for i in 1:length(p_out)]..., G...]
            if criterion == 1
                is_out = create_is_out_1(qcp, sub_interval)
            end
            if criterion == 2
                f_bounds = bounds(qcp.problem.f, X_0, sub_interval)
                is_out = create_is_out_2(qcp, sub_interval, f_bounds, 0.00001)
            end
            is_out_union |= is_out(X_0)
            increment!(indices, lengths, indices_forall)
            if isempty(indices_forall)
                break
            end
        end
        for i in indices_forall
            indices[i] = 1
        end
        is_out_intersection &= is_out_union
        increment!(indices, lengths, indices_exists)
        if isempty(indices_exists)
            break
        end
    end
    return is_out_intersection
end

check_is_out_1(X_0, p_in, G, qcp) = check_is_out(X_0, p_in, G, qcp, 1)
check_is_out_2(X_0, p_in, G, qcp) = check_is_out(X_0, p_in, G, qcp, 2)

function pave(X::IntervalArithmetic.IntervalBox{N, T}, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, check_is_in, check_is_out)::Tuple{Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}} where {N, T<:Number}
    inn = []
    @assert nand(allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) "Refinement and subdivision are mutually exclusive. Use --help for more information."
    p_in_0 = deepcopy(p_in)
    p_out_0 = deepcopy(p_out)
    inn = []
    out = []
    delta = []
    list = [(X, p_in, p_out)]
    while !isempty(list)
        X, p_in, p_out = pop!(list)
        if !allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection
            if check_is_in(X, p_in, G, qcp)
                push!(inn, X)
            elseif check_is_out(X, p_out, G, qcp)
                push!(out, X)
            elseif IntervalArithmetic.diam(X) < ϵ_x
                push!(delta, X)
            else
                X_1, X_2 = bisect(X)
                push!(list, (X_1, deepcopy(p_in_0), deepcopy(p_out_0)))
                push!(list, (X_2, deepcopy(p_in_0), deepcopy(p_out_0)))
            end
        end
        if allow_exists_and_forall_bisection
            p_in_max = maximum(IntervalArithmetic.diam.(first.(p_in)); init=0.0)
            p_out_max = maximum(IntervalArithmetic.diam.(first.(p_out)); init=0.0)
            p_max = maximum((p_in_max, p_out_max))
            X_max = IntervalArithmetic.diam(X)
            if check_is_in(X, p_in, G, qcp)
                push!(inn, X)
            elseif check_is_out(X, p_out, G, qcp)
                push!(out, X)
            elseif X_max < ϵ_x && p_max < ϵ_p
                push!(delta, X)
            elseif X_max >= ϵ_x && p_max < ϵ_p
                X_1, X_2 = bisect(X)
                push!(list, (X_1, deepcopy(p_in_0), deepcopy(p_out_0)))
                push!(list, (X_2, deepcopy(p_in_0), deepcopy(p_out_0)))
            elseif X_max < ϵ_x && p_max >= ϵ_p
                if ϵ_p <= p_in_max
                    bisect_largest_forall!(p_in, qcp.qvs, qcp.p, qcp.n)
                end
                if ϵ_p <= p_out_max
                    bisect_largest_exists!(p_out, qcp.qvs, qcp.p, qcp.n)
                end
                push!(list, (X, p_in, p_out))
            else
                if X_max < p_max
                    if ϵ_p <= p_in_max
                        bisect_largest_forall!(p_in, qcp.qvs, qcp.p, qcp.n)
                    end
                    if ϵ_p <= p_out_max
                        bisect_largest_exists!(p_out, qcp.qvs, qcp.p, qcp.n)
                    end
                    push!(list, (X, p_in, p_out))
                else
                    X_1, X_2 = bisect(X)
                    push!(list, (X_1, deepcopy(p_in_0), deepcopy(p_out_0)))
                    push!(list, (X_2, deepcopy(p_in_0), deepcopy(p_out_0)))
                end
            end
        end
        if allow_exists_or_forall_bisection
            indices_forall = [i for (q, i) in qcp.qvs if q == Forall] .- length(X)
            indices_exists = [i for (q, i) in qcp.qvs if q == Exists] .- length(X)
            p_in_max = maximum(IntervalArithmetic.diam.(first.(p_in[indices_exists])); init=0.0)
            p_out_max = maximum(IntervalArithmetic.diam.(first.(p_out[indices_forall])); init=0.0)
            p_max = maximum((p_in_max, p_out_max))
            X_max = IntervalArithmetic.diam(X)
            if check_is_in(X, p_in, G, qcp)
                push!(inn, X)
            elseif check_is_out(X, p_out, G, qcp)
                push!(out, X)
            elseif X_max < ϵ_x && p_max < ϵ_p
                push!(delta, X)
            elseif X_max >= ϵ_x && p_max < ϵ_p
                X_1, X_2 = bisect(X)
                push!(list, (X_1, deepcopy(p_in_0), deepcopy(p_out_0)))
                push!(list, (X_2, deepcopy(p_in_0), deepcopy(p_out_0)))
            elseif X_max < ϵ_x && p_max >= ϵ_p
                if ϵ_p <= p_in_max
                    bisect_largest_forall!(p_in, qcp.qvs, qcp.p, qcp.n)
                end
                if ϵ_p <= p_out_max
                    bisect_largest_exists!(p_out, qcp.qvs, qcp.p, qcp.n)
                end
                push!(list, (X, p_in, p_out))
            else
                if X_max < p_max
                    if ϵ_p <= p_in_max
                        bisect_largest_exists!(p_in, qcp.qvs, qcp.p, qcp.n)
                    end
                    if ϵ_p <= p_out_max
                        bisect_largest_forall!(p_out, qcp.qvs, qcp.p, qcp.n)
                    end
                    push!(list, (X, p_in, p_out))
                else
                    X_1, X_2 = bisect(X)
                    push!(list, (X_1, deepcopy(p_in_0), deepcopy(p_out_0)))
                    push!(list, (X_2, deepcopy(p_in_0), deepcopy(p_out_0)))
                end
            end
        end
    end
    return inn, out, delta
end

pave_11(X, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) = pave(X, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, check_is_in_1, check_is_out_1)
pave_12(X, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) = pave(X, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, check_is_in_1, check_is_out_2)
pave_21(X, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) = pave(X, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, check_is_in_2, check_is_out_1)
pave_22(X, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) = pave(X, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, check_is_in_2, check_is_out_2)

# Utils

function volume_box(box)
    return prod(IntervalArithmetic.diam.(box))
end

function volume_boxes(boxes)
    if isempty(boxes)
        return 0
    end
    return sum(volume_box.(boxes))
end

function merge_intervals(intervals)
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

## Plots drawing

rectangle(p, q) = Shape([p[1],q[1],q[1],p[1]], [p[2],p[2],q[2],q[2]])

function draw_lines(pl, intervals, color)
    for interval in intervals
        p = [interval.lo, -0.1]
        q = [interval.hi, 0.1]
        plot!(pl, rectangle(p,q), color=color, linecolor=nothing, legend=:false)
    end
end

function draw_rows(p, boxes, color)
    draw_lines(p, merge_intervals([box[1] for box in boxes]), color)
end

draw_inn_lines(p, inns) = draw_rows(p, inns, :green)
draw_out_lines(p, outs) = draw_rows(p, outs, :cyan)
draw_delta_lines(p, deltas) = draw_rows(p, deltas, :yellow)

function draw_rectangles(boxes, color)
    for box in boxes
        x_interval = box[1]
        y_interval = box[2]
        p = (x_interval.lo, y_interval.lo)
        q = (x_interval.hi, y_interval.hi)
        plot!(rectangle(p, q), color=color, legend=:false)
    end
end

draw_delta_rectangles(delta) = draw_rectangles(delta, :yellow)
draw_inn_rectangles(inn) = draw_rectangles(inn, :green)
draw_out_rectangles(out) = draw_rectangles(out, :cyan)

function draw(p, X_0, inn, out, delta)
    if isa(X_0, IntervalBox{1, <:Number})
        xticks!((X_0[1].lo:1:X_0[1].hi))
        yaxis!(false)

        draw_delta_lines(p, delta)
        draw_inn_lines(p, inn)
        draw_out_lines(p, out)
    elseif isa(X_0, IntervalBox{2, <:Number})
        xs = X_0[1]
        ys = X_0[2]
        xlims!((xs.lo, xs.hi))
        ylims!((ys.lo, ys.hi))

        draw_delta_rectangles(delta)
        draw_inn_rectangles(inn)
        draw_out_rectangles(out)
    else
        error("Plotting is only supported for 1D and 2D problems.")
    end
end

function print_inn_out_delta(inn, out, delta)
    Base.println("Union of elements of inn: ", merge_intervals([box[1] for box in inn]))
    Base.println("Union of elements of out: ", merge_intervals([box[1] for box in out]))
    Base.println("Union of elements of delta: ", merge_intervals([box[1] for box in delta]))
end

function print_delta_width(delta)
    Base.println("Width of delta regions: ", IntervalArithmetic.diam.(merge_intervals([box[1] for box in delta])))
end

## Luxor drawing

function luxor_box2pq(box)
    x = box[1]
    y = box[2]
    p = Point(x.lo, y.lo)
    q = Point(x.hi, y.hi)
    return p, q
end

function luxor_rescale(p, q, X_0, width, height, buffer)
    if isa(X_0, IntervalBox{1, <:Number})
        scale_x = width / (X_0[1].hi - X_0[1].lo)
        scale_y = height / 0.2
        p_rescaled = Point(buffer + (p.x - X_0[1].lo) * scale_x, buffer + height - (p.y + 0.1) * scale_y)
        q_rescaled = Point(buffer + (q.x - X_0[1].lo) * scale_x, buffer + height - (q.y + 0.1) * scale_y)
    elseif isa(X_0, IntervalBox{2, <:Number})
        scale_x = width / (X_0[1].hi - X_0[1].lo)
        scale_y = height / (X_0[2].hi - X_0[2].lo)
        p_rescaled = Point(buffer + (p.x - X_0[1].lo) * scale_x, buffer + height - (p.y - X_0[2].lo) * scale_y)
        q_rescaled = Point(buffer + (q.x - X_0[1].lo) * scale_x, buffer + height - (q.y - X_0[2].lo) * scale_y)
    end
    return p_rescaled, q_rescaled
end

function luxor_rescaled_pq(box, X_0, width, height, buffer)
    p, q = luxor_box2pq(box)
    return luxor_rescale(p, q, X_0, width, height, buffer)
end

function luxor_draw_rows(boxes, color, X_0, width, height, buffer)
    sethue(color)
    for box in boxes
        p = Point(box[1].lo, -0.1)
        q = Point(box[1].hi, 0.1)
        p_rescaled, q_rescaled = luxor_rescale(p, q, X_0, width, height, buffer)
        Luxor.box(p_rescaled, q_rescaled, :fill)
    end
end

luxor_draw_inn_rows(inn, X_0, width, height, buffer) = luxor_draw_rows(inn, "green", X_0, width, height, buffer)
luxor_draw_out_rows(out, X_0, width, height, buffer) = luxor_draw_rows(out, "cyan", X_0, width, height, buffer)
luxor_draw_delta_rows(delta, X_0, width, height, buffer) = luxor_draw_rows(delta, "yellow", X_0, width, height, buffer)

function luxor_draw_boxes(boxes, color, X_0, width, height, buffer)
    sethue(color)
    for box in boxes
        p, q = luxor_rescaled_pq(box, X_0, width, height, buffer)
        Luxor.box(p, q, :fill)
    end
end

luxor_draw_inn_boxes(inn, X_0, width, height, buffer) = luxor_draw_boxes(inn, "green", X_0, width, height, buffer)
luxor_draw_out_boxes(out, X_0, width, height, buffer) = luxor_draw_boxes(out, "cyan", X_0, width, height, buffer)
luxor_draw_delta_boxes(delta, X_0, width, height, buffer) = luxor_draw_boxes(delta, "yellow", X_0, width, height, buffer)

function luxor_draw(X_0, inn, out, delta, width, height, buffer)
    if isa(X_0, IntervalBox{1, <:Number})
        background("white")

        luxor_draw_inn_rows(inn, X_0, width, height, buffer)
        luxor_draw_out_rows(out, X_0, width, height, buffer)
        luxor_draw_delta_rows(delta, X_0, width, height, buffer)

        sethue("black")
        # xticks
        tickline(Point(buffer, buffer + height), Point(buffer + width, buffer + height), startnumber= X_0[1].lo, finishnumber=X_0[1].hi, major=4, minor=0)
    elseif isa(X_0, IntervalBox{2, <:Number})
        background("white")

        luxor_draw_inn_boxes(inn, X_0, width, height, buffer)
        luxor_draw_out_boxes(out, X_0, width, height, buffer)
        luxor_draw_delta_boxes(delta, X_0, width, height, buffer)

        sethue("black")
        # xticks
        tickline(Point(buffer, buffer + height), Point(buffer + width, buffer + height), startnumber= X_0[1].lo, finishnumber=X_0[1].hi, major=4, minor=0)
        # yticks
        tickline(Point(buffer + width, buffer + height), Point(buffer + width, buffer), startnumber= X_0[2].lo, finishnumber=X_0[2].hi, major=4, minor=0)
    else
        error("Plotting is only supported for 1D and 2D problems.")
    end
end