using Plots
using Luxor
using MathTeXEngine

include("genreach2.jl")
include("quantifiedconstraintproblem.jl")

const global minus_inf = -100000
const global plus_inf = 100000
const global strict_epsilon = 0.0001

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
    return bisect(box, i)..., i
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

function complement(x::IntervalArithmetic.Interval{T}) where T <: Number
    l = []
    if x.lo != minus_inf
        push!(l, interval(minus_inf, x.lo - strict_epsilon))
    end
    if x.hi != plus_inf
        push!(l, interval(x.hi + strict_epsilon, plus_inf))
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

"Custom version of isbounded to handle pseudo_infinity values."
function is_bounded(interval::IntervalArithmetic.Interval{T}) where T <: Number
    return !isempty(interval) && interval.lo != minus_inf && interval.hi != plus_inf
end

function is_bounded(intervals::Vector{IntervalArithmetic.Interval{T}}) where T <: Number
    return all(is_bounded, intervals)
end

function atom_conjunction(atom::Atom, n::Int)
    if n == 1
        return atom
    end
    conjuncted_predicates = [Atom("$(atom.value)_$(i)") for i in 1:n]
    return ∧(conjuncted_predicates...)
end

function atom_complement_disjunction(atom::Atom, intervals::Vector{IntervalArithmetic.Interval{T}}) where T <: Number
    if length(intervals) == 1
        interval = first(intervals)
        if is_bounded(interval)
            atom_1 = Atom("$(atom.value)_1")
            atom_2 = Atom("$(atom.value)_2")
            return atom_1 ∨ atom_2
        else
            return atom
        end
    end

    i = 1
    disjuncted_predicates = []
    for interval in intervals
        value_i = "$(atom.value)_$(i)"
        if is_bounded(interval)
            atom_1 = Atom(value_i * "_1")
            atom_2 = Atom(value_i * "_2")
            push!(disjuncted_predicates, atom_1 ∨ atom_2)
        else
            push!(disjuncted_predicates, Atom(value_i))
        end
        i += 1
    end
    if length(disjuncted_predicates) == 1
        return first(disjuncted_predicates)
    else
        return ∨(disjuncted_predicates...)
    end
end

function expand_1(M, indices_dict, intervals::Vector{IntervalArithmetic.Interval{T}}, atoms) where T <: Number
    if M isa Atom
        indices = indices_dict[M.value]
        n_indices = last(indices) - first(indices) + 1
        return atom_conjunction(M, n_indices)
    elseif token(M) == ∧
        return expand_1(first(M.children), indices_dict, intervals, atoms) ∧ expand_1(last(M.children), indices_dict, intervals, atoms)
    elseif token(M) == ∨
        return expand_1(first(M.children), indices_dict, intervals, atoms) ∨ expand_1(last(M.children), indices_dict, intervals, atoms)
    elseif token(M) == ¬
        atom = first(M.children)
        if !(atom isa Atom)
            error("The formula should only contain negations at the leaves.")
        else
            range = indices_dict[atom.value]
            sub_intervals = intervals[range]
            return atom_complement_disjunction(atom, sub_intervals)
        end
    else
        error("Unknown syntax tree type: $(typeof(M)).")
    end
end

expand_1(M, indices_dict, intervals) = expand_1(M, indices_dict, intervals, atoms(M))

function expand_2(M, indices_dict, intervals::Vector{IntervalArithmetic.Interval{T}}, atoms, previous_token=nothing) where T <: Number
    if M isa Atom
        if previous_token != ¬
            range = indices_dict[M.value]
            sub_intervals = intervals[range]
            return atom_complement_disjunction(M, sub_intervals)
        else
            indices = indices_dict[M.value]
            n_indices = last(indices) - first(indices) + 1
            return atom_conjunction(M, n_indices)
        end
    elseif token(M) == ∧
        if previous_token == ¬
            error("The formula should only contain negations at the leaves.")
        else
            return expand_2(first(M.children), indices_dict, intervals, atoms) ∨ expand_2(last(M.children), indices_dict, intervals, atoms)
        end
    elseif token(M) == ∨
        if previous_token == ¬
            error("The formula should only contain negations at the leaves.")
        else
            return expand_2(first(M.children), indices_dict, intervals, atoms) ∧ expand_2(last(M.children), indices_dict, intervals, atoms)
        end
    elseif token(M) == ¬
        return expand_2(first(M.children), indices_dict, intervals, atoms, ¬)
    else
        error("Unknown syntax tree type: $(typeof(M)).")
    end
end

expand_2(M, indices_dict, intervals) = expand_2(M, indices_dict, intervals, atoms(M))

# function expand_negations_1(M, indices_dict, intervals, atoms)
#     if M isa Atom
#         return M
#     elseif token(M) == ∧
#         return expand_negations_1(first(M.children), intervals, atoms) ∧ expand_negations_1(last(M.children), intervals, atoms)
#     elseif token(M) == ∨
#         return expand_negations_1(first(M.children), intervals, atoms) ∨ expand_negations_1(last(M.children), intervals, atoms)
#     elseif token(M) == ¬
#         atom = first(M.children)
#         pos = findfirst(x -> x == atom, atoms)
#         interval = intervals[pos]
#         if is_bounded(interval)
#             atom_1 = Atom(atom.value * "_1")
#             atom_2 = Atom(atom.value * "_2")
#             return atom_1 ∨ atom_2
#         else
#             return atom
#         end
#     else
#         error("Unknown syntax tree type: $(typeof(M))")
#     end
# end

# expand_negations_1(M, intervals) = expand_negations_1(M, intervals, atoms(M))

# function expand_negations_2(M, intervals, indices_dict, atoms, previous_token=nothing)
#     if M isa Atom
#         if previous_token != ¬
#             pos = findfirst(x -> x == M, atoms)
#             interval = intervals[pos]
#             if is_bounded(interval)
#                 atom_1 = Atom(M.value * "_1")
#                 atom_2 = Atom(M.value * "_2")
#                 return atom_1 ∨ atom_2
#             else
#                 return M
#             end
#         else
#             return M
#         end
#     elseif token(M) == ∧
#         return expand_negations_2(first(M.children), intervals, atoms) ∨ expand_negations_2(last(M.children), intervals, atoms)
#     elseif token(M) == ∨
#         return expand_negations_2(first(M.children), intervals, atoms) ∧ expand_negations_2(last(M.children), intervals, atoms)
#     elseif token(M) == ¬
#         return expand_negations_2(first(M.children), intervals, atoms, ¬)
#     else
#         error("Unknown syntax tree type: $(typeof(M))")
#     end
# end

# expand_negations_2(M, intervals) = expand_negations_2(M, intervals, atoms(M))

function duplication_positions_1(M, indices_dict, intervals)
    positions = Int[]
    to_visit = Union{Atom,SyntaxBranch}[M]
    while !isempty(to_visit)
        current = pop!(to_visit)
        if current isa Atom
            continue
        elseif (token(current) == ∧) || (token(current) == ∨)
            append!(to_visit, children(current))
        elseif token(current) == ¬
            atom = first(current.children)
            range = indices_dict[atom.value]
            for i in range
                interval = intervals[i]
                if is_bounded(interval)
                    push!(positions, i)
                end
            end
        else
            error("Unknown syntax tree type: $(typeof(current))")
        end
    end
    return positions
end

function duplication_positions_2(M, indices_dict, intervals)
    positions = Int[]
    to_visit = Tuple{Union{Atom,SyntaxBranch},Union{Nothing,Connective}}[(M, nothing)]
    while !isempty(to_visit)
        current, previous_token = pop!(to_visit)
        if current isa Atom
            if previous_token != ¬
                range = indices_dict[current.value]
                for i in range
                    interval = intervals[i]
                    if is_bounded(interval)
                        push!(positions, i)
                    end
                end
            end
        elseif (token(current) == ∧) || (token(current) == ∨)
            child_1, child_2 = children(current)
            push!(to_visit, (child_1, nothing))
            push!(to_visit, (child_2, nothing))
        elseif token(current) == ¬
        else
            error("Unknown syntax tree type: $(typeof(current))")
        end
    end
    return positions
end

function complement_ranges_1(M, indices_dict)
    ranges = UnitRange[]
    to_visit = Union{Atom,SyntaxBranch}[M]
    while !isempty(to_visit)
        current = pop!(to_visit)
        if current isa Atom
            continue
        elseif (token(current) == ∧) || (token(current) == ∨)
            append!(to_visit, children(current))
        elseif token(current) == ¬
            atom = first(current.children)
            range = indices_dict[atom.value]
            push!(ranges, range)
        else
            error("Unknown syntax tree type: $(typeof(current))")
        end
    end
    return ranges
end

function complement_ranges_2(M, indices_dict)
    ranges = UnitRange[]
    to_visit = Tuple{Union{Atom,SyntaxBranch},Union{Nothing,Connective}}[(M, nothing)]
    while !isempty(to_visit)
        current, previous_token = pop!(to_visit)
        if current isa Atom
            if previous_token != ¬
                range = indices_dict[current.value]
                push!(ranges, range)
            end
        elseif (token(current) == ∧) || (token(current) == ∨)
            child_1, child_2 = children(current)
            push!(to_visit, (child_1, nothing))
            push!(to_visit, (child_2, nothing))
        elseif token(current) == ¬
            continue
        else
            error("Unknown syntax tree type: $(typeof(current))")
        end
    end
    return ranges
end

function expand_intervals(intervals, ranges)
    positions = vcat(map(collect, ranges)...)
    expanded_intervals = []
    for i in 1:length(intervals)
        if i in positions
            complement_intervals = complement(intervals[i])
            for interval in complement_intervals
                push!(expanded_intervals, interval)
            end
        else
            push!(expanded_intervals, intervals[i])
        end
    end
    return expanded_intervals
end

function inflate(l, indices)
    inflated = []
    for (i, item) in enumerate(l)
        if i in indices
            push!(inflated, item)
            push!(inflated, item)
        else
            push!(inflated, item)
        end
    end
    return inflated
end

function is_zero_in(R, i)
    return !isempty(R[i]) && interval(0, 0) ⊆ interval(min(R[i]), max(R[i]))
end

function is_zero_not_in(R, i)
    return isempty(R[i]) || interval(0,0) ⊈ interval(min(R[i]), max(R[i]))
end

function test_zero_in(M, R, positions_dict)
    @match M begin
        ::Atom => is_zero_in(R, positions_dict[M.value])
        ::SyntaxBranch where (token(M) == ∧) => all(test_zero_in(child, R, positions_dict) for child in M.children)
        ::SyntaxBranch where (token(M) == ∨) => any(test_zero_in(child, R, positions_dict) for child in M.children)
        ::SyntaxBranch where (token(M) == ¬) => error("The formula should not contain negations at this point. Negations should have been expanded.")
        _ => error("Unknown syntax tree type: $(typeof(M))")
    end
end

function test_zero_not_in(M, R, positions_dict)
    @match M begin
        ::Atom => is_zero_not_in(R, positions_dict[M.value])
        ::SyntaxBranch where (token(M) == ∧) => any(test_zero_not_in(child, R, positions_dict) for child in M.children)
        ::SyntaxBranch where (token(M) == ∨) => all(test_zero_not_in(child, R, positions_dict) for child in M.children)
        ::SyntaxBranch where (token(M) == ¬) => error("The formula should not contain negations at this point. Negations should have been expanded.")
        _ => error("Unknown syntax tree type: $(typeof(M))")
    end
end

function create_is_in_1(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Forall, i) for i in 1:length(X)]..., qcp.qvs..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        quantifiers_relaxed = [[[(Forall, i) for i in 1:length(X)]..., qcp.qvs_relaxed[j]..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...] for j in 1:qcp.n]
        dirty_qs = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        problem = qcp.problem
        # G = last(intervals, qcp.n)
        R_inner = QEapprox_o0_inner(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qcp.p, qcp.n, [X.v..., intervals...])
        return test_zero_in(problem.M, R_inner, problem.positions_dict)
    end
end

function create_is_in_2(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Exists, i) for i in 1:length(X)]..., negation.(qcp.qvs)..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        quantifiers_relaxed = [[[(Exists, i) for i in 1:length(X)]..., negation.(qcp.qvs_relaxed[j])..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...] for j in 1:qcp.n]
        dirty_qs = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        problem = qcp.problem
        # G_minus = [interval(-∞, intervals[end-i].lo) ∩ f_bounds[end-i] for i in (qcp.n-1):-1:0]
        # G = last(intervals, qcp.n)
        R_outer = QEapprox_o0_outer(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qcp.p, qcp.n, [X.v..., intervals...])
        return test_zero_not_in(problem.M, R_outer, problem.positions_dict) 
    end
end

function create_is_out_1(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Exists, i) for i in 1:length(X)]..., qcp.qvs..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        quantifiers_relaxed = [[[(Exists, i) for i in 1:length(X)]..., qcp.qvs_relaxed[j]..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...] for j in 1:qcp.n]
        dirty_qs = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        problem = qcp.problem
        R_outer = QEapprox_o0_outer(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qcp.p, qcp.n, [X.v..., intervals...])
        return test_zero_not_in(problem.M, R_outer, problem.positions_dict)
    end
end

function create_is_out_2(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where {N, T<:Number}
        quantifiers = [[(Forall, i) for i in 1:length(X)]..., negation.(qcp.qvs)..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        quantifiers_relaxed = [[[(Forall, i) for i in 1:length(X)]..., negation.(qcp.qvs_relaxed[j])..., [(Exists, qcp.p-i) for i in (qcp.n-1):-1:0]...] for j in 1:qcp.n]
        dirty_qs = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        problem = qcp.problem
        # G_complement = complement_disjunction(last(intervals, qcp.n), problem.dnf_indices)
        # G = last(intervals, qcp.n)
        R_inner = QEapprox_o0_inner(problem.f, problem.Df, dirty_quantifiers, dirty_qs, qcp.p, qcp.n, [X.v..., intervals...])
        return test_zero_in(problem.M, R_inner, problem.positions_dict)
    end
end

# function bounds(f, X, interval)
#     return [f[i]([X.v..., interval...]) for i in 1:length(f)]
# end

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
            if criterion == 1
                is_out = create_is_out_1(qcp, sub_interval)
            end
            if criterion == 2
                is_out = create_is_out_2(qcp, sub_interval)
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

global z_in, z_out

function pave(X::IntervalArithmetic.IntervalBox{N, T}, p_in, p_out, G, variables, sizes, f_num, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, criterion_in, criterion_out)::Tuple{Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}} where {N, T<:Number}
    @assert ((allow_exists_and_forall_bisection || allow_exists_or_forall_bisection) && !isnothing(ϵ_p)) || (!allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection) "ϵ_p must be provided when bisection on parameter space is allowed."
    @assert nand(allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) "Refinement and subdivision are mutually exclusive. Use --help for more information."
    @assert length(G) == qcp.n "Length of G must be equal to the number of functions, n = $(qcp.n)."
    @assert length(X) + length(p_in) == qcp.p "Total number of variables, in X and p_in, must be equal to p = $(qcp.p)."
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
    check_is_in = criterion_in == 1 ? check_is_in_1 : check_is_in_2
    check_is_out = criterion_out == 1 ? check_is_out_1 : check_is_out_2
    expand_in = criterion_in == 1 ? expand_1 : expand_2
    expand_out = criterion_out == 1 ? expand_1 : expand_2
    # expand_negations_in = criterion_in == 1 ? expand_negations_1 : expand_negations_2
    # expand_negations_out = criterion_out == 1 ? expand_negations_1 : expand_negations_2
    duplication_positions_in = criterion_in == 1 ? duplication_positions_1 : duplication_positions_2
    duplication_positions_out = criterion_out == 1 ? duplication_positions_1 : duplication_positions_2
    complement_ranges_in = criterion_in == 1 ? complement_ranges_1 : complement_ranges_2
    complement_ranges_out = criterion_out == 1 ? complement_ranges_1 : complement_ranges_2
    indices_dict = get_indices_dict(qcp.problem.M, sizes)
    M_in = expand_in(qcp.problem.M, indices_dict, G)
    M_out = expand_out(qcp.problem.M, indices_dict, G)
    # M_out = parseformula("P_1_1 ∧ P_1_2")
    # display(M_in)
    # display(M_out)
    # exit()
    positions_in = duplication_positions_in(qcp.problem.M, indices_dict, G)
    positions_out = duplication_positions_out(qcp.problem.M, indices_dict, G)
    n = qcp.n
    n_in = length(positions_in)
    n_out = length(positions_out)
    @variables z_in[n + n_in]
    @variables z_out[n + n_out]
    f_num_in = inflate(f_num, positions_in)
    for i in 1:(n + n_in)
        f_num_in[i] -= z_in[i]
    end
    f_num_out = inflate(f_num, positions_out)
    for i in 1:(n + n_out)
        f_num_out[i] -= z_out[i]
    end
    f_fun_in, Df_fun_in = build_function_f_Df(f_num_in, [variables..., z_in...], n + n_in, qcp.p + n + n_in)
    f_fun_out, Df_fun_out = build_function_f_Df(f_num_out, [variables..., z_out...], n + n_out, qcp.p + n + n_out)
    G_in = expand_intervals(G, complement_ranges_in(qcp.problem.M, indices_dict))
    G_out = expand_intervals(G, complement_ranges_out(qcp.problem.M, indices_dict))
    # display(G_in)
    # display(G_out)
    # exit()
    qvs_relaxed_in = inflate(qcp.qvs_relaxed, positions_in)
    qvs_relaxed_out = inflate(qcp.qvs_relaxed, positions_out)
    sizes_in = inflate(sizes, positions_in)
    sizes_out = inflate(sizes, positions_out)
    positions_dict_in = get_positions_dict(M_in)
    positions_dict_out = get_positions_dict(M_out)
    problem_in = Problem(M_in, f_fun_in, Df_fun_in, positions_dict_in)
    problem_out = Problem(M_out, f_fun_out, Df_fun_out, positions_dict_out)
    qcp_in = QuantifiedConstraintProblem(problem_in, qcp.qvs, qvs_relaxed_in, qcp.p + n + n_in, n + n_in)
    qcp_out = QuantifiedConstraintProblem(problem_out, qcp.qvs, qvs_relaxed_out, qcp.p + n + n_out, n + n_out)
    p_in_0 = deepcopy(p_in)
    p_out_0 = deepcopy(p_out)
    inn = []
    out = []
    delta = []
    list = [(X, p_in, p_out)]
    while !isempty(list)
        X, p_in, p_out = pop!(list)
        if !allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection
            if check_is_in(X, p_in, G_in, qcp_in)
                push!(inn, X)
            elseif check_is_out(X, p_out, G_out, qcp_out)
                push!(out, X)
            elseif all(map(<, IntervalArithmetic.diam.(X), ϵ_x))
                push!(delta, X)
            else
                X_1, X_2, _ = bisect_precision(X, ϵ_x)
                push!(list, (X_1, deepcopy(p_in_0), deepcopy(p_out_0)))
                push!(list, (X_2, deepcopy(p_in_0), deepcopy(p_out_0)))
            end
        end
        if allow_exists_and_forall_bisection
            p_in_diams = IntervalArithmetic.diam.(first.(p_in))
            p_out_diams = IntervalArithmetic.diam.(first.(p_out))
            p_maxs = max.(p_in_diams, p_out_diams)
            X_diams = IntervalArithmetic.diam.(X)
            if check_is_in(X, p_in, G_in, qcp_in)
                push!(inn, X)
            elseif check_is_out(X, p_out, G_out, qcp_out)
                push!(out, X)
            elseif all(map(<, X_diams, ϵ_x)) && all(map(<, p_maxs, ϵ_p))
                push!(delta, X)
            elseif any(map(>=, X_diams, ϵ_x)) && all(map(<, p_maxs, ϵ_p))
                X_1, X_2, _ = bisect_precision(X, ϵ_x)
                push!(list, (X_1, deepcopy(p_in_0), deepcopy(p_out_0)))
                push!(list, (X_2, deepcopy(p_in_0), deepcopy(p_out_0)))
            elseif all(map(<, X_diams, ϵ_x)) && any(map(>=, p_maxs, ϵ_p))
                if any(map(<=, ϵ_p, p_in_diams))
                    bisect_largest_forall!(p_in, qcp.qvs, qcp.p, qcp.n, ϵ_p)
                    # bisect_eps_forall!(p_in, qcp.qvs, ϵ_p, qcp.p, qcp.n)
                end
                if any(map(<=, ϵ_p, p_out_diams))
                    bisect_largest_exists!(p_out, qcp.qvs, qcp.p, qcp.n, ϵ_p)
                    # bisect_eps_exists!(p_out, qcp.qvs, ϵ_p, qcp.p, qcp.n)
                end
                push!(list, (X, p_in, p_out))
            else
                if maximum(X_diams) < maximum(p_maxs)
                    if any(map(<=, ϵ_p, p_in_diams))
                        bisect_largest_forall!(p_in, qcp.qvs, qcp.p, qcp.n, ϵ_p)
                        # bisect_eps_forall!(p_in, qcp.qvs, ϵ_p, qcp.p, qcp.n)
                    end
                    if any(map(<=, ϵ_p, p_out_diams))
                        bisect_largest_exists!(p_out, qcp.qvs, qcp.p, qcp.n, ϵ_p)
                        # bisect_eps_exists!(p_out, qcp.qvs, ϵ_p, qcp.p, qcp.n)
                    end
                    push!(list, (X, p_in, p_out))
                else
                    X_1, X_2, _ = bisect_precision(X, ϵ_x)
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
            if check_is_in(X, p_in, G_in, qcp_in)
                push!(inn, X)
            elseif check_is_out(X, p_out, G_out, qcp_out)
                push!(out, X)
            elseif all(map(<, X_diams, ϵ_x)) && all(map(<, p_maxs, ϵ_p))
                push!(delta, X)
            elseif any(map(>=, X_diams, ϵ_x)) && all(map(<, p_maxs, ϵ_p))
                X_1, X_2, _ = bisect_precision(X, ϵ_x)
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
                    X_1, X_2, _ = bisect_precision(X, ϵ_x)
                    push!(list, (X_1, deepcopy(p_in_0), deepcopy(p_out_0)))
                    push!(list, (X_2, deepcopy(p_in_0), deepcopy(p_out_0)))
                end
            end
        end
    end
    return inn, out, delta
end

pave_11(X, p_in, p_out, G, variables, sizes, f_num, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) = pave(X, p_in, p_out, G, variables, sizes, f_num, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, 1, 1)
pave_12(X, p_in, p_out, G, variables, sizes, f_num, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) = pave(X, p_in, p_out, G, variables, sizes, f_num, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, 1, 2)
pave_21(X, p_in, p_out, G, variables, sizes, f_num, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) = pave(X, p_in, p_out, G, variables, sizes, f_num, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, 2, 1)
pave_22(X, p_in, p_out, G, variables, sizes, f_num, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) = pave(X, p_in, p_out, G, variables, sizes, f_num, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, 2, 2)

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

function bisect_increasing_order(box, i)
    box_1, box_2 = bisect(box, i)
    return box_1, box_2
end

function bisect_decreasing_order(box, i)
    box_1, box_2 = bisect(box, i)
    return box_2, box_1
end

function create_bisect_in_optimization_direction(sign)
    if sign > 0
        return bisect_increasing_order
    else
        return bisect_decreasing_order
    end
end

"""
Pave the input box X, when each component of X is monotonous with respect to the set membership relation.
Faster than pave_monotonous_sides, but relies on a heuristic that may produce additional undecided boxes. (This behavior is obvious when the input domain X_0 is totally inside or outside the set.)
"""
function pave_monotonous_mid(X::IntervalArithmetic.IntervalBox{N, T}, optimization_directions, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, check_is_in, check_is_out) where {N, T<:Number}
    @assert ((allow_exists_and_forall_bisection || allow_exists_or_forall_bisection) && !isnothing(ϵ_p)) || (!allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection) "ϵ_p must be provided when bisection on parameter space is allowed."
    @assert nand(allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) "Refinement and subdivision are mutually exclusive. Use --help for more information."
    @assert length(optimization_directions) == length(X) "Length of optimization_directions must be equal to the number of variables in X."
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
    queue = [(X, p_in, p_out)]
    while !isempty(queue)
        X, p_in, p_out = pop!(queue)
        if !allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection
            error("TO DO")
        end
        if allow_exists_and_forall_bisection
            if all(map(<, IntervalArithmetic.diam.(X), ϵ_x))
                push!(delta, X)
                continue
            end
            dim_slice, X_slice = bisection_slice(X, ϵ_x)
            if check_is_in(X_slice, p_in, G, qcp)
                bisect_in_optimization_direction = create_bisect_in_optimization_direction(optimization_directions[dim_slice])
                X_inn, X_queue = bisect_in_optimization_direction(X, dim_slice)
                push!(inn, X_inn)
                push!(queue, (X_queue, p_in, p_out))
            elseif check_is_out(X_slice, p_out, G, qcp)
                bisect_in_optimization_direction = create_bisect_in_optimization_direction(optimization_directions[dim_slice])
                X_queue, X_out = bisect_in_optimization_direction(X, dim_slice)
                push!(out, X_out)
                push!(queue, (X_queue, p_in, p_out))
            else
                X_1, X_2, _ = bisect_precision(X, ϵ_x)
                push!(queue, (X_1, p_in, p_out))
                push!(queue, (X_2, p_in, p_out))
            end
        end
        if allow_exists_or_forall_bisection
            error("TO DO")
        end
    end
    return inn, out, delta
end

function backward(interval, sign)
    if sign > 0
        return interval.lo
    else
        return interval.hi
    end
end

function forward(interval, sign)
    if sign > 0
        return interval.hi
    else
        return interval.lo
    end
end

function slice(box, optimization_directions, dim, directed_function)
    l = []

    for i in 1:length(box)
        if i == dim
            directed_value = directed_function(box[dim], optimization_directions[dim])
            push!(l, interval(directed_value, directed_value))
        else
            push!(l, box[i])
        end
    end

    return IntervalBox(l)
end

backward_slice(box, optimization_directions, dim) = slice(box, optimization_directions, dim, backward)
forward_slice(box, optimization_directions, dim) = slice(box, optimization_directions, dim, forward)

function slices(box, optimization_directions, unchanged_face, directed_slicing_function)
    l = []
    for i in 1:length(box)
        if !isnothing(unchanged_face) && i == dim(unchanged_face) && bound(directed_slicing_function, optimization_directions[i]) == bound(unchanged_face)
            continue
        end
        push!(l, directed_slicing_function(box, optimization_directions, i))
    end
    return l
end

@enum Bound begin
    Upper
    Lower
end

Face = Tuple{Int, Bound}

function dim(face::Face)
    return face[1]
end

function bound(face::Face)
    return face[2]
end

function bound(f, sign)
    if f == forward_slice && sign == 1
        return Upper
    end
    if f == backward_slice && sign == -1
        return Upper
    end
    if f == forward_slice && sign == -1
        return Lower
    end
    if f == backward_slice && sign == 1
        return Lower
    end
end

backward_slices(box, optimization_directions, unchanged_face) = slices(box, optimization_directions, unchanged_face, backward_slice)
forward_slices(box, optimization_directions, unchanged_face) = slices(box, optimization_directions, unchanged_face, forward_slice)

"""
Pave the input box X, when each component of X is monotonous with respect to the set membership relation.
Slower than pave_monotonous_mid, but does not produce unexpected undecided boxes.
"""
function pave_monotonous_sides_optimized(X::IntervalArithmetic.IntervalBox{N, T}, optimization_directions, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, check_is_in, check_is_out) where {N, T<:Number}
    @assert ((allow_exists_and_forall_bisection || allow_exists_or_forall_bisection) && !isnothing(ϵ_p)) || (!allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection) "ϵ_p must be provided when bisection on parameter space is allowed."
    @assert nand(allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) "Refinement and subdivision are mutually exclusive. Use --help for more information."
    @assert length(optimization_directions) == length(X) "Length of optimization_directions must be equal to the number of variables in X."
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
    stack = Tuple{IntervalArithmetic.IntervalBox, Union{Face, Nothing}, Vector{Vector{IntervalArithmetic.Interval}}, Vector{Vector{IntervalArithmetic.Interval}}}[]
    push!(stack, (X, nothing, p_in, p_out))
    while !isempty(stack)
        X, unchanged_face, p_in, p_out = popfirst!(stack)
        if !allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection
            error("TO DO")
        end
        if allow_exists_and_forall_bisection
            backward_X_slices = backward_slices(X, optimization_directions, unchanged_face)
            forward_X_slices = forward_slices(X, optimization_directions, unchanged_face)
            if any(check_is_in(X_slice, p_in, G, qcp) for X_slice in forward_X_slices)
                push!(inn, X)
            elseif any(check_is_out(X_slice, p_out, G, qcp) for X_slice in backward_X_slices)
                push!(out, X)
            elseif all(map(<, IntervalArithmetic.diam.(X), ϵ_x))
                push!(delta, X)
            else
                X_1, X_2, dim = bisect_precision(X, ϵ_x)
                unchanged_face_1 = (dim, Lower)
                unchanged_face_2 = (dim, Upper)
                push!(stack, (X_1, unchanged_face_1, p_in, p_out))
                push!(stack, (X_2, unchanged_face_2, p_in, p_out))
            end
        end
        if allow_exists_or_forall_bisection
            error("TO DO")
        end
    end
    return inn, out, delta
end

pave_monotonous_mid_12(X, optimization_directions, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) = pave_monotonous_mid(X, optimization_directions, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, check_is_in_1, check_is_out_2)
pave_monotonous_sides_optimized_12(X, optimization_directions, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) = pave_monotonous_sides_optimized(X, optimization_directions, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection, check_is_in_1, check_is_out_2)

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
    p = Luxor.Point(x.lo, y.lo)
    q = Luxor.Point(x.hi, y.hi)
    return p, q
end

function luxor_rescale(p, q, X_0, width, height, buffer)
    if isa(X_0, IntervalBox{1, <:Number})
        scale_x = width / (X_0[1].hi - X_0[1].lo)
        scale_y = height / 0.2
        p_rescaled = Luxor.Point(buffer + (p.x - X_0[1].lo) * scale_x, buffer + height - (p.y + 0.1) * scale_y)
        q_rescaled = Luxor.Point(buffer + (q.x - X_0[1].lo) * scale_x, buffer + height - (q.y + 0.1) * scale_y)
    elseif isa(X_0, IntervalBox{2, <:Number})
        scale_x = width / (X_0[1].hi - X_0[1].lo)
        scale_y = height / (X_0[2].hi - X_0[2].lo)
        p_rescaled = Luxor.Point(buffer + (p.x - X_0[1].lo) * scale_x, buffer + height - (p.y - X_0[2].lo) * scale_y)
        q_rescaled = Luxor.Point(buffer + (q.x - X_0[1].lo) * scale_x, buffer + height - (q.y - X_0[2].lo) * scale_y)
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
        p = Luxor.Point(box[1].lo, -0.1)
        q = Luxor.Point(box[1].hi, 0.1)
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
        tickline(Luxor.Point(buffer, buffer + height), Luxor.Point(buffer + width, buffer + height), startnumber= X_0[1].lo, finishnumber=X_0[1].hi, major=4, minor=0)
    elseif isa(X_0, IntervalBox{2, <:Number})
        background("white")

        luxor_draw_inn_boxes(inn, X_0, width, height, buffer)
        luxor_draw_out_boxes(out, X_0, width, height, buffer)
        luxor_draw_delta_boxes(delta, X_0, width, height, buffer)

        sethue("black")
        # xticks
        tickline(Luxor.Point(buffer, buffer + height), Luxor.Point(buffer + width, buffer + height), startnumber= X_0[1].lo, finishnumber=X_0[1].hi, major=4, minor=0)
        # yticks
        tickline(Luxor.Point(buffer + width, buffer + height), Luxor.Point(buffer + width, buffer), startnumber= X_0[2].lo, finishnumber=X_0[2].hi, major=4, minor=0)
    else
        error("Plotting is only supported for 1D and 2D problems.")
    end
end