# using IntervalArithmetic
using Plots
using TimerOutputs

const to = TimerOutput()

include("genreach2.jl")
include("quantifiedconstraintproblem.jl")

const global minus_inf = -100000
const global plus_inf = 100000

# Paving

function bisect_eps(interval, ϵ)
    parts = [interval]
    while diam(first(parts)) >= ϵ
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
        intervals[i] = bisect_eps(intervals[i][1], eps[i])
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

function bisect_largest_quantifier!(intervals, qvs, p, n, quantifier, ϵ)
    pos_quantifier = [i for (q, i) in qvs if q == quantifier] .- (p - n - length(qvs))
    diams = [if (i in pos_quantifier) IntervalArithmetic.diam(first(intervals[i])) else -1.0 end for i in 1:length(intervals)]
    is_not_bisectable = diams .< ϵ
    diams[is_not_bisectable] .= -1.0
    pos_max = argmax(diams)
    parts = []
    for interval in intervals[pos_max]
        a, b = IntervalArithmetic.bisect(interval)
        push!(parts, a)
        push!(parts, b)
    end
    intervals[pos_max] = parts
end

bisect_largest_exists!(intervals, qvs, p, n, ϵ) = bisect_largest_quantifier!(intervals, qvs, p, n, Exists, ϵ)
bisect_largest_forall!(intervals, qvs, p, n, ϵ) = bisect_largest_quantifier!(intervals, qvs, p, n, Forall, ϵ)

function bisect_precision(box, ϵ)
    diams = IntervalArithmetic.diam.(box)
    is_not_bisectable = diams .< ϵ
    copy_diams = [d for d in diams]
    copy_diams[is_not_bisectable] .= -1.0
    i = argmax(copy_diams)
    return bisect(box, i)
end

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

function complement(x::IntervalArithmetic.Interval; ϵ=0.0001)
    l = []
    if x.lo != minus_inf
        push!(l, interval(minus_inf, x.lo - ϵ))
    end
    if x.hi != plus_inf
        push!(l, interval(x.hi + ϵ, plus_inf))
    end
    return l
end

struct Component
    interval::IntervalArithmetic.Interval
    index::Int
end

function complement_conjunction_components(G::Vector{IntervalArithmetic.Interval{T}}, conjunction_indices::Vector{Int}) where T <: Number
    l = []
    for (i, index) in enumerate(conjunction_indices)
        compl_G_i = complement(G[index])
        for sub_compl_G_i in compl_G_i
            sub_region = [Component(interval(minus_inf, plus_inf), j) for j in conjunction_indices]
            sub_region[i] = Component(sub_compl_G_i, index)
            push!(l, sub_region)
        end
    end
    return l
end

function complement_disjunction(G::Vector{IntervalArithmetic.Interval{T}}, dnf_indices::Vector{Vector{Int}}) where T <: Number
    conjunctions = []
    for conjunction_indices in dnf_indices
        compl = complement_conjunction_components(G, conjunction_indices)
        push!(conjunctions, compl)
    end

    n = length(G)
    conjunction_product = collect(Iterators.product(conjunctions...))
    disjunction = []
    for components in conjunction_product
        comp = reduce(vcat, components)
        filtered_intervals = []
        for i in 1:n
            ith_intervals = [c.interval for c in comp if c.index == i]
            ith_interval =  reduce(intersect, ith_intervals; init=interval(minus_inf, plus_inf))
            push!(filtered_intervals, ith_interval)
        end
        push!(disjunction, filtered_intervals)
    end
    return disjunction
end

function disjunction(G::Vector{IntervalArithmetic.Interval{T}}, dnf_indices::Vector{Vector{Int}}) where T <: Number
    n = length(G)
    disjunction = []
    for conjunction_indices in dnf_indices
        sub_region = repeat([interval(minus_inf, plus_inf)], n)
        for i in conjunction_indices
            sub_region[i] = G[i]
        end
        push!(disjunction, sub_region)
    end
    return disjunction
end

function create_is_in_1(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Forall, i) for i in 1:length(X)]..., qcp.qvs..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        quantifiers_relaxed = [[[(Forall, i) for i in 1:length(X)]..., qcp.qvs_relaxed[j]..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...] for j in 1:qcp.n]
        dirty_qs = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        problem = qcp.problem
        G = disjunction(last(intervals, qcp.n), problem.dnf_indices)
        for G_i in G
            @timeit to "in" R_inner, _ = QEapprox_o0(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qcp.p, qcp.n, [X.v..., intervals[1:end-qcp.n]..., G_i...])
            if all(!isempty(R_inner[i]) && interval(0, 0) ⊆ interval(min(R_inner[i]), max(R_inner[i])) for i in 1:qcp.n)
                return true
            end
        end
        return false
    end
end

function create_is_in_2(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Exists, i) for i in 1:length(X)]..., negation.(qcp.qvs)..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        quantifiers_relaxed = [[[(Exists, i) for i in 1:length(X)]..., negation.(qcp.qvs_relaxed[j])..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...] for j in 1:qcp.n]
        dirty_qs = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        problem = qcp.problem
        G_minus = [interval(-∞, intervals[end-i].lo) ∩ ϕ_bounds[end-i] for i in (qcp.n-1):-1:0]
        if any(isempty, G_minus)
            test_minus = true
        else
            _, R_outer_minus = QEapprox_o0(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qcp.p, qcp.n, [X.v..., intervals[1:end-qcp.n]..., G_minus...])
            test_minus = any([interval(0, 0) ⊈ interval(min(R_outer_minus[i]), max(R_outer_minus[i])) for i in 1:qcp.n])
        end
        G_plus = [interval(intervals[end-i].hi, ∞) ∩ ϕ_bounds[end-i] for i in (qcp.n-1):-1:0]
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
        # G = [intervals[end-i] ∩ ϕ_bounds[end-i] for i in (qcp.n-1):-1:0]
        # if any(isempty, G)
        #     return true
        # end
        # _, R_outer = QEapprox_o0(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qcp.p, qcp.n, [X.v..., intervals[1:end-qcp.n]..., G...])
        _, R_outer = QEapprox_o0(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qcp.p, qcp.n, [X.v..., intervals...])
        if all(isempty, R_outer)
            return true
        end
        return  any([interval(0,0) ⊈ interval(min(R_outer[i]), max(R_outer[i])) for i in 1:qcp.n if !isempty(R_outer[i])])
    end
end

function create_is_out_2(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}}, ϵ::Float64=0.0001)::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Forall, i) for i in 1:length(X)]..., negation.(qcp.qvs)..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        quantifiers_relaxed = [[[(Forall, i) for i in 1:length(X)]..., negation.(qcp.qvs_relaxed[j])..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...] for j in 1:qcp.n]
        dirty_qs = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        problem = qcp.problem
        G_complement = complement_disjunction(last(intervals, qcp.n), problem.dnf_indices)
        for G_complement_i in G_complement
            @timeit to "out" R_inner, _ = QEapprox_o0(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qcp.p, qcp.n, [X.v..., intervals[1:end-qcp.n]..., G_complement_i...])
            if all(!isempty(R_inner[i]) && interval(0, 0) ⊆ interval(min(R_inner[i]), max(R_inner[i])) for i in 1:qcp.n)
                return true
            end
        end
        return false
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
                is_in = create_is_in_2(qcp, sub_interval)
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
            # f_bounds = bounds(qcp.problem.f, X_0, sub_interval)
            if criterion == 1
                is_out = create_is_out_1(qcp, sub_interval)
            end
            if criterion == 2
                is_out = create_is_out_2(qcp, sub_interval, 0.00001)
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
    @assert ((allow_exists_and_forall_bisection || allow_exists_or_forall_bisection) && !isnothing(ϵ_p)) || (!allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection) "ϵ_p must be provided when bisection on parameter space is allowed."
    @assert nand(allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) "Refinement and subdivision are mutually exclusive. Use --help for more information."
    @assert length(G) == qcp.n "Length of G must be equal to the number of functions, n = $(qcp.n)."
    @assert length(X) + length(p_in) + length(G) == qcp.p "Total number of variables, in X and p_in, must be equal to p - n = $(qcp.p - qcp.n)."
    @assert length(qcp.qvs) == length(p_in) "Number of quantified variables must be equal to the number of parameter boxes, $(length(p_in))."
    for qv in qcp.qvs
        @assert length(X) < index(qv) <= length(X) + length(p_in) "Quantified variables must be in the parameter space: indices between $(length(X)+1) and $(qcp.p - qcp.n)."
    end
    for qvs in qcp.qvs_relaxed
        @assert length(qvs) == length(p_in) "Number of quantified variables must be equal to the number of parameter boxes, $(length(p_in))."
        for qv in qvs
            @assert length(X) < index(qv) <= length(X) + length(p_in) "Quantified variables must be in the parameter space: indices between $(length(X)+1) and $(qcp.p - qcp.n)."
        end
    end
    inn = []
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
            elseif all(map(<, IntervalArithmetic.diam.(X), ϵ_x))
                push!(delta, X)
            else
                X_1, X_2 = bisect_precision(X, ϵ_x)
                push!(list, (X_1, deepcopy(p_in_0), deepcopy(p_out_0)))
                push!(list, (X_2, deepcopy(p_in_0), deepcopy(p_out_0)))
            end
        end
        if allow_exists_and_forall_bisection
            p_in_diams = IntervalArithmetic.diam.(first.(p_in))
            p_out_diams = IntervalArithmetic.diam.(first.(p_out))
            p_maxs = max.(p_in_diams, p_out_diams)
            X_diams = IntervalArithmetic.diam.(X)
            if check_is_in(X, p_in, G, qcp)
                push!(inn, X)
            elseif check_is_out(X, p_out, G, qcp)
                push!(out, X)
            elseif all(map(<, X_diams, ϵ_x)) && all(map(<, p_maxs, ϵ_p))
                push!(delta, X)
            elseif any(map(>=, X_diams, ϵ_x)) && all(map(<, p_maxs, ϵ_p))
                X_1, X_2 = bisect_precision(X, ϵ_x)
                push!(list, (X_1, deepcopy(p_in_0), deepcopy(p_out_0)))
                push!(list, (X_2, deepcopy(p_in_0), deepcopy(p_out_0)))
            elseif all(map(<, X_diams, ϵ_x)) && any(map(>=, p_maxs, ϵ_p))
                if any(map(<=, ϵ_p, p_in_diams))
                    # bisect_largest_forall!(p_in, qcp.qvs, qcp.p, qcp.n, ϵ_p)
                    bisect_eps_forall!(p_in, qcp.qvs, ϵ_p, qcp.p, qcp.n)
                end
                if any(map(<=, ϵ_p, p_out_diams))
                    # bisect_largest_exists!(p_out, qcp.qvs, qcp.p, qcp.n, ϵ_p)
                    bisect_eps_exists!(p_out, qcp.qvs, ϵ_p, qcp.p, qcp.n)
                end
                push!(list, (X, p_in, p_out))
            else
                if maximum(X_diams) < maximum(p_maxs)
                    if any(map(<=, ϵ_p, p_in_diams))
                        # bisect_largest_forall!(p_in, qcp.qvs, qcp.p, qcp.n, ϵ_p)
                        bisect_eps_forall!(p_in, qcp.qvs, ϵ_p, qcp.p, qcp.n)
                    end
                    if any(map(<=, ϵ_p, p_out_diams))
                        # bisect_largest_exists!(p_out, qcp.qvs, qcp.p, qcp.n, ϵ_p)
                        bisect_eps_exists!(p_out, qcp.qvs, ϵ_p, qcp.p, qcp.n)
                    end
                    push!(list, (X, p_in, p_out))
                else
                    X_1, X_2 = bisect_precision(X, ϵ_x)
                    push!(list, (X_1, deepcopy(p_in_0), deepcopy(p_out_0)))
                    push!(list, (X_2, deepcopy(p_in_0), deepcopy(p_out_0)))
                end
            end
        end
        if allow_exists_or_forall_bisection
            indices_forall = [i for (q, i) in qcp.qvs if q == Forall] .- length(X)
            indices_exists = [i for (q, i) in qcp.qvs if q == Exists] .- length(X)
            p_in_diams = IntervalArithmetic.diam.(first.(p_in))
            p_in_diams[indices_forall] .= -1.0
            p_out_diams = IntervalArithmetic.diam.(first.(p_out))
            p_out_diams[indices_exists] .= -1.0
            p_maxs = max.(p_in_diams, p_out_diams)
            X_diams = IntervalArithmetic.diam.(X)
            if check_is_in(X, p_in, G, qcp)
                push!(inn, X)
            elseif check_is_out(X, p_out, G, qcp)
                push!(out, X)
            elseif all(map(<, X_diams, ϵ_x)) && all(map(<, p_maxs, ϵ_p))
                push!(delta, X)
            elseif any(map(>=, X_diams, ϵ_x)) && all(map(<, p_maxs, ϵ_p))
                X_1, X_2 = bisect_precision(X, ϵ_x)
                push!(list, (X_1, deepcopy(p_in_0), deepcopy(p_out_0)))
                push!(list, (X_2, deepcopy(p_in_0), deepcopy(p_out_0)))
            elseif all(map(<, X_diams, ϵ_x)) && any(map(>=, p_maxs, ϵ_p))
                if any(map(<=, ϵ_p, p_in_diams))
                    bisect_largest_forall!(p_in, qcp.qvs, qcp.p, qcp.n, ϵ_p)
                end
                if any(map(<=, ϵ_p, p_out_diams))
                    bisect_largest_exists!(p_out, qcp.qvs, qcp.p, qcp.n, ϵ_p)
                end
                push!(list, (X, p_in, p_out))
            else
                if maximum(X_diams) < maximum(p_maxs)
                    if any(map(<=, ϵ_p, p_in_diams))
                        bisect_largest_exists!(p_in, qcp.qvs, qcp.p, qcp.n, ϵ_p)
                    end
                    if any(map(<=, ϵ_p, p_out_diams))
                        bisect_largest_forall!(p_out, qcp.qvs, qcp.p, qcp.n, ϵ_p)
                    end
                    push!(list, (X, p_in, p_out))
                else
                    X_1, X_2 = bisect_precision(X, ϵ_x)
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

function bisection_slice(box, ϵ)
    diams = IntervalArithmetic.diam.(box)
    is_not_bisectable = diams .< ϵ
    copy_diams = [d for d in diams]
    copy_diams[is_not_bisectable] .= -1.0
    pos_max = argmax(copy_diams)
    mid_max = mid(box[pos_max])
    l = []
    for i in 1:length(box)
        if i != pos_max
            push!(l, box[i])
        else
            push!(l, interval(mid_max, mid_max))
        end
    end
    return pos_max, IntervalBox(l)
end

function decrease(box, i)
    box_1, box_2 = bisect(box, i)
    return box_1
end

function increase(box, i)
    box_1, box_2 = bisect(box, i)
    return box_2
end

function monotony(signs, i)
    if signs[i] > 0
        return increase
    else
        return decrease
    end
end

function antimonotony(signs, i)
    if signs[i] > 0
        return decrease
    else
        return increase
    end
end

function pave_monotonous(X::IntervalArithmetic.IntervalBox{N, T}, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, check_is_in, check_is_out) where {N, T<:Number}
    @assert ((allow_exists_and_forall_bisection || allow_exists_or_forall_bisection) && !isnothing(ϵ_p)) || (!allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection) "ϵ_p must be provided when bisection on parameter space is allowed."
    @assert nand(allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) "Refinement and subdivision are mutually exclusive. Use --help for more information."
    @assert length(G) == qcp.n "Length of G must be equal to the number of functions, n = $(qcp.n)."
    @assert length(X) + length(p_in) + length(G) == qcp.p "Total number of variables, in X and p_in, must be equal to p - n = $(qcp.p - qcp.n)."
    @assert length(qcp.qvs) == length(p_in) "Number of quantified variables must be equal to the number of parameter boxes, $(length(p_in))."
    for qv in qcp.qvs
        @assert length(X) < index(qv) <= length(X) + length(p_in) "Quantified variables must be in the parameter space: indices between $(length(X)+1) and $(qcp.p - qcp.n)."
    end
    for qvs in qcp.qvs_relaxed
        @assert length(qvs) == length(p_in) "Number of quantified variables must be equal to the number of parameter boxes, $(length(p_in))."
        for qv in qvs
            @assert length(X) < index(qv) <= length(X) + length(p_in) "Quantified variables must be in the parameter space: indices between $(length(X)+1) and $(qcp.p - qcp.n)."
        end
    end
    inn = []
    p_in_0 = deepcopy(p_in)
    p_out_0 = deepcopy(p_out)
    p_in_diams= IntervalArithmetic.diam.(first.(p_in))
    p_out_diams = IntervalArithmetic.diam.(first.(p_out))
    if any(map(<=, ϵ_p, p_in_diams))
        # bisect_largest_forall!(p_in, qcp.qvs, qcp.p, qcp.n)
        bisect_eps_forall!(p_in, qcp.qvs, ϵ_p, qcp.p, qcp.n)
    end
    if any(map(<=, ϵ_p, p_out_diams))
        # bisect_largest_exists!(p_out, qcp.qvs, qcp.p, qcp.n)
        bisect_eps_exists!(p_out, qcp.qvs, ϵ_p, qcp.p, qcp.n)
    end
    inn = []
    out = []
    delta = []
    X_in = X
    X_out = deepcopy(X)
    list = [(X, p_in, p_out)]
    while !isempty(list)
        X, p_in, p_out = pop!(list)
        if !allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection
            error("TO DO")
        end
        if allow_exists_and_forall_bisection
            if all(map(<, IntervalArithmetic.diam.(X), ϵ_x))
                push!(delta, X)
                continue
            end
            X_mid = IntervalBox(mid(X))
            if check_is_in(X_mid, p_in, G, qcp)
                X_1, X_2 = bisect_precision(X, ϵ_x)
                push!(inn, X_2)
                push!(list, (X_1, p_in, p_out))
            elseif check_is_out(X_mid, p_out, G, qcp)
                X_1, X_2 = bisect_precision(X, ϵ_x)
                push!(out, X_1)
                push!(list, (X_2, p_in, p_out))
            else
                X_1, X_2 = bisect_precision(X, ϵ_x)
                push!(list, (X_1, p_in, p_out))
                push!(list, (X_2, p_in, p_out))
            end
        end
        if allow_exists_or_forall_bisection
            error("TO DO")
        end
    end
    return inn, out, delta
end

pave_monotonous_12(X, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) = pave_monotonous(X, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, check_is_in_1, check_is_out_2)

function pave_monotonous_2D(X::IntervalArithmetic.IntervalBox{N, T}, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, check_is_in, check_is_out) where {N, T<:Number}
    @assert ((allow_exists_and_forall_bisection || allow_exists_or_forall_bisection) && !isnothing(ϵ_p)) || (!allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection) "ϵ_p must be provided when bisection on parameter space is allowed."
    @assert nand(allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) "Refinement and subdivision are mutually exclusive. Use --help for more information."
    @assert length(G) == qcp.n "Length of G must be equal to the number of functions, n = $(qcp.n)."
    @assert length(X) + length(p_in) + length(G) == qcp.p "Total number of variables, in X and p_in, must be equal to p - n = $(qcp.p - qcp.n)."
    @assert length(qcp.qvs) == length(p_in) "Number of quantified variables must be equal to the number of parameter boxes, $(length(p_in))."
    for qv in qcp.qvs
        @assert length(X) < index(qv) <= length(X) + length(p_in) "Quantified variables must be in the parameter space: indices between $(length(X)+1) and $(qcp.p - qcp.n)."
    end
    for qvs in qcp.qvs_relaxed
        @assert length(qvs) == length(p_in) "Number of quantified variables must be equal to the number of parameter boxes, $(length(p_in))."
        for qv in qvs
            @assert length(X) < index(qv) <= length(X) + length(p_in) "Quantified variables must be in the parameter space: indices between $(length(X)+1) and $(qcp.p - qcp.n)."
        end
    end
    inn = []
    p_in_0 = deepcopy(p_in)
    p_out_0 = deepcopy(p_out)
    p_in_max = maximum(IntervalArithmetic.diam.(first.(p_in)); init=0.0)
    p_out_max = maximum(IntervalArithmetic.diam.(first.(p_out)); init=0.0)
    if ϵ_p <= p_in_max
        # bisect_largest_forall!(p_in, qcp.qvs, qcp.p, qcp.n)
        bisect_eps_forall!(p_in, qcp.qvs, ϵ_p, qcp.p, qcp.n)
    end
    if ϵ_p <= p_out_max
        # bisect_largest_exists!(p_out, qcp.qvs, qcp.p, qcp.n)
        bisect_eps_exists!(p_out, qcp.qvs, ϵ_p, qcp.p, qcp.n)
    end
    inn = []
    out = []
    delta = []
    X_in = X
    X_out = deepcopy(X)
    list = [(X_in, X_out, p_in, p_out)]
    while !isempty(list)
        X_in, X_out, p_in, p_out = pop!(list)
        if !allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection
            error("TO DO")
        end
        if allow_exists_and_forall_bisection
            if !isempty(X_in)
                X_in_mid = IntervalBox(mid(X_in))
                if check_is_in(X_in_mid, p_in, G, qcp)
                    push!(inn, X_in_mid[1])
                    if IntervalArithmetic.diam(X_in) >= ϵ_x
                        new_X_in, _ = bisect(X_in)
                    else
                        new_X_in = IntervalBox(emptyinterval())
                    end
                elseif IntervalArithmetic.diam(X_in) >= ϵ_x
                    _, new_X_in = bisect(X_in)
                else
                    new_X_in = IntervalBox(emptyinterval())
                end
            else
                new_X_in = IntervalBox(emptyinterval())
            end
            if !isempty(X_out)
                X_out_mid = IntervalBox(mid(X_out))
                if check_is_out(X_out_mid, p_out, G, qcp)
                    push!(out, X_out_mid[1])
                    if IntervalArithmetic.diam(X_out) >= ϵ_x
                        _, new_X_out = bisect(X_out)
                    else
                        new_X_out = IntervalBox(emptyinterval())
                    end
                elseif IntervalArithmetic.diam(X_out) >= ϵ_x
                    new_X_out, _ = bisect(X_out)
                else
                    new_X_out = IntervalBox(emptyinterval())
                end
            else
                new_X_out = IntervalBox(emptyinterval())
            end
            if !isempty(new_X_in) || !isempty(new_X_out)
                push!(list, (new_X_in, new_X_out, p_in, p_out))
            end
        end
        if allow_exists_or_forall_bisection
            error("TO DO")
        end
    end
    return inn, out, delta
end

pave_monotonous_2D_12(X, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) = pave_monotonous_2D(X, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, check_is_in_1, check_is_out_2)

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