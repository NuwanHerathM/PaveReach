using IntervalArithmetic
# using TimerOutputs
using Plots

# const to = TimerOutput()

include("genreach2.jl")
include("quantifierproblem.jl")

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
    @assert sum(length.(intervals)) ==  length(intervals) "Each interval should be a single interval."
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

function bisect_largest!(intervals)
    (_, pos_max) = findmax(IntervalArithmetic.diam.(first.(intervals)))
    parts = []
    for interval in intervals[pos_max]
        a, b = IntervalArithmetic.bisect(interval)
        push!(parts, a)
        push!(parts, b)
    end
    intervals[pos_max] = parts
end

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

function create_is_in_1(qe::QuantifierProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Forall, i) for i in 1:length(X)]..., qe.qvs..., [(Exists, qe.p) for i in 1:qe.n]...]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        qs = [[[(Forall, i) for i in 1:length(X)]..., qe.q[j]..., [(Exists, qe.p) for i in 1:qe.n]...] for j in 1:qe.n]
        dirty_qs = quantifiedvariables2dirtyvariables.(qs)
        problem = qe.problem
        R_inner, _ = QEapprox_o0(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qe.p, qe.n, [X.v..., intervals...])
        # println(R_inner)
        if any(isempty, R_inner)
            return false
        end
        return all(interval(0, 0) ⊆ interval(min(R_inner[i]), max(R_inner[i])) for i in 1:qe.n)
    end
end

function create_is_out_1(qe::QuantifierProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Exists, i) for i in 1:length(X)]..., qe.qvs..., [(Exists, qe.p) for i in 1:qe.n]...]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        qs = [[[(Exists, i) for i in 1:length(X)]..., qe.q[j]..., [(Exists, qe.p) for i in 1:qe.n]...] for j in 1:qe.n]
        dirty_qs = quantifiedvariables2dirtyvariables.(qs)
        problem = qe.problem
        _, R_outer = QEapprox_o0(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qe.p, qe.n, [X.v..., intervals...])
        # println(R_outer)
        if all(isempty, R_outer)
            return true
        end
        return  any([interval(0,0) ⊈ interval(min(R_outer[i]), max(R_outer[i])) for i in 1:qe.n if !isempty(R_outer[i])])
    end
end

function check_is_in(X_0, p_in, G, qe)
    pg_in = [p_in..., G...]

    indices_forall = [i for (q, i) in qe.qvs if q == Forall] .- length(X_0)
    indices_exists = [i for (q, i) in qe.qvs if q == Exists] .- length(X_0)
    append!(indices_exists, (length(p_in)+1):(length(p_in)+length(G)))

    indices = [1 for i in 1:(length(p_in)+length(G))]
    lengths_pg_in = length.(pg_in)
    # println(lengths_pg_in)
    is_in_union = false
    while !is_in_union && (isempty(indices_exists) || indices[indices_exists] <= lengths_pg_in[indices_exists])
        is_in_intersection = true
        while is_in_intersection && (isempty(indices_forall) || indices[indices_forall] <= lengths_pg_in[indices_forall])
            sub_interval = [pg_in[i][indices[i]] for i in 1:length(pg_in)]
            is_in = create_is_in_1(qe, sub_interval)
            is_in_intersection &= is_in(X_0)
            if isempty(indices_forall)
                break
            end
            increment!(indices, lengths_pg_in, indices_forall)
        end
        for i in indices_forall
            indices[i] = 1
        end
        is_in_union |= is_in_intersection
        if isempty(indices_exists)
            break
        end
        increment!(indices, lengths_pg_in, indices_exists)
    end
    return is_in_union
end


function check_is_out(X_0, p_out, G, qe)
    pg_out = [p_out..., G...]

    indices_forall = [i for (q, i) in qe.qvs if q == Forall] .- length(X_0)
    indices_exists = [i for (q, i) in qe.qvs if q == Exists] .- length(X_0)
    append!(indices_exists, (length(p_in)+1):(length(p_in)+length(G)))

    indices = [1 for i in 1:(length(p_in)+length(G))]
    lengths_pg_out = length.(pg_out)
    # println(lengths_pg_out)
    is_out_intersection = true
    while is_out_intersection && (isempty(indices_exists) || indices[indices_exists] <= lengths_pg_out[indices_exists])
        is_out_union = false
        while !is_out_union && (isempty(indices_forall) || indices[indices_forall] <= lengths_pg_out[indices_forall])
            sub_interval = [pg_out[i][indices[i]] for i in 1:length(pg_out)]
            is_out = create_is_out_1(qe, sub_interval)
            is_out_union |= is_out(X_0)
            increment!(indices, lengths_pg_out, indices_forall)
            if isempty(indices_forall)
                break
            end
        end
        for i in indices_forall
            indices[i] = 1
        end
        is_out_intersection &= is_out_union
        increment!(indices, lengths_pg_out, indices_exists)
        if isempty(indices_exists)
            break
        end
    end
    return is_out_intersection
end

function pave(X, p_in, p_out, G, qe, ϵ, is_refined, allow_normal_p_bisect)
    @assert nand(is_refined, allow_normal_p_bisect) "Cannot have P refined and normal bisection on P."
    p_in_0 = deepcopy(p_in)
    p_out_0 = deepcopy(p_out)
    G_0 = deepcopy(G)
    inn = []
    out = []
    delta = []
    list = [(X, p_in, p_out, G)]
    while !isempty(list)
        X, p_in, p_out, G = pop!(list)
        if check_is_in(X, p_in, G, qe)
            push!(inn, X)
        elseif check_is_out(X, p_out, G, qe)
            push!(out, X)
        elseif IntervalArithmetic.diam(X) < ϵ
            push!(delta, X)
        else
            if is_refined
                p_in_max = maximum(IntervalArithmetic.diam.(first.(p_in)))
                p_out_max = maximum(IntervalArithmetic.diam.(first.(p_out)))
                p_max = maximum((p_in_max, p_out_max))
                G_max = maximum(IntervalArithmetic.diam.(first.(G)))
                X_max = IntervalArithmetic.diam(X)
                if X_max < p_max && p_max < G_max
                    bisect_largest!(G)
                    push!(list, (X, p_in, p_out, G))
                    continue
                end
                if X_max < p_max
                    # bisect_largest!(p_in)
                    # bisect_largest!(p_out)
                    bisect_largest_forall!(p_in, qe.qvs, qe.p, qe.n)
                    bisect_largest_exists!(p_out, qe.qvs, qe.p, qe.n)
                    push!(list, (X, p_in, p_out, G))
                    continue
                end
            end
            if allow_normal_p_bisect
                indices_forall = [i for (q, i) in qe.qvs if q == Forall] .- length(X)
                indices_exists = [i for (q, i) in qe.qvs if q == Exists] .- length(X)
                p_in_max = isempty(indices_exists) ? -1 : maximum(IntervalArithmetic.diam.(first.(p_in[indices_exists])))
                p_out_max = isempty(indices_forall) ? -1 : maximum(IntervalArithmetic.diam.(first.(p_out[indices_forall])))
                p_max = maximum((p_in_max, p_out_max))
                G_max = maximum(IntervalArithmetic.diam.(first.(G)))
                X_max = IntervalArithmetic.diam(X)
                if X_max < p_max && p_max < G_max
                    bisect_largest!(G)
                    push!(list, (X, p_in, p_out, G))
                    continue
                end
                if X_max < p_max
                    bisect_largest_exists!(p_in, qe.qvs, qe.p, qe.n)
                    bisect_largest_forall!(p_out, qe.qvs, qe.p, qe.n)
                    push!(list, (X, p_in, p_out, G))
                    continue
                end
            end
            X_1, X_2 = bisect(X)
            push!(list, (X_1, deepcopy(p_in_0), deepcopy(p_out_0), deepcopy(G_0)))
            push!(list, (X_2, deepcopy(p_in_0), deepcopy(p_out_0), deepcopy(G_0)))
        end
    end
    return inn, out, delta
end

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

rectangle(p, q) = Shape([p[1],q[1],q[1],p[1]], [p[2],p[2],q[2],q[2]])

function draw_lines(intervals, color, y)
    for interval in intervals
        p = [interval.lo, y-0.1]
        q = [interval.hi, y+0.1]
        plot!(rectangle(p,q), color=color, linecolor=nothing, legend=:false)
    end
end

function draw_rows(blocks, color)
    n = length(blocks)
    for i in 1:n
        y = n - i
        draw_lines(merge_intervals([box[1] for box in blocks[i]]), color, y)
    end
end

draw_inn_lines(inns) = draw_rows(inns, :green)
draw_out_lines(outs) = draw_rows(outs, :cyan)
draw_delta_lines(deltas) = draw_rows(deltas, :yellow)

function draw(X_0, inn, out, delta)
    default(size=(600,60), grid=false)
    xlabel!("x")
    xticks!((X_0[1].lo:1:X_0[1].hi))
    yaxis!(false)

    draw_delta_lines(delta)
    draw_inn_lines(inn)
    draw_out_lines(out)
end