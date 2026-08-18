using Plots
using Luxor
using MathTeXEngine
using LinearAlgebra

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

function bisect_eps_quantifier!(intervals, qvs, eps, p, quantifier)
    @assert sum(length.(intervals); init=0) ==  length(intervals) "Each interval should be a single interval."
    pos_quantifier = [i for (q, i) in qvs if q == quantifier] .- (p - length(qvs))

    for i in pos_quantifier
        intervals[i] = bisect_eps(intervals[i][1], eps[i])
    end
end

bisect_eps_exists!(intervals, qvs, eps, p) = bisect_eps_quantifier!(intervals, qvs, eps, p, Exists)
bisect_eps_forall!(intervals, qvs, eps, p) = bisect_eps_quantifier!(intervals, qvs, eps, p, Forall)

function pointify_quantifier!(intervals, qvs, p, quantifier)
    pos_quantifier = [i for (q, i) in qvs if q == quantifier] .- (p - length(qvs))

    for i in pos_quantifier
        intervals[i] = interval.(mid.(intervals[i]))
    end
end

pointify_exists!(intervals, qvs, p) = pointify_quantifier!(intervals, qvs, p, Exists)
pointify_forall!(intervals, qvs, p) = pointify_quantifier!(intervals, qvs, p, Forall)

function refine_in!(P_in, qvs, eps, p)
    bisect_eps_exists!(P_in, qvs, eps, p)
    pointify_exists!(P_in, qvs, p)
end

function refine_out!(P_out, qvs, eps, p)
    bisect_eps_forall!(P_out, qvs, eps, p)
    pointify_forall!(P_out, qvs, p)
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

function bisect_largest_quantifier!(intervals, qvs, p, quantifier, ϵ)
    pos_quantifier = [i for (q, i) in qvs if q == quantifier] .- (p - length(qvs))
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

bisect_largest_exists!(intervals, qvs, p, ϵ) = bisect_largest_quantifier!(intervals, qvs, p, Exists, ϵ)
bisect_largest_forall!(intervals, qvs, p, ϵ) = bisect_largest_quantifier!(intervals, qvs, p, Forall, ϵ)

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
    ranges = UnitRange{Int}[]
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
    ranges = UnitRange{Int}[]
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

function expand_intervals(intervals::Vector{IntervalArithmetic.Interval{T}}, ranges::Vector{UnitRange{Int}}) where T <: Number
    positions = vcat(map(collect, ranges)...)
    expanded_intervals = IntervalArithmetic.Interval{T}[]
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
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where N
        quantifiers = [[(Forall, i) for i in 1:length(X)]; get_qvs(qcp); [(Exists, get_p(qcp)-i) for i in (get_n(qcp)-1):-1:0]]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        qvs_relaxed = get_qvs_relaxed(qcp)
        quantifiers_relaxed = [[[(Forall, i) for i in 1:length(X)]; qvs_relaxed[j]; [(Exists, get_p(qcp)-i) for i in (get_n(qcp)-1):-1:0]] for j in 1:get_n(qcp)]
        dirty_qs = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        R_inner = QEapprox_o0_inner(get_f(qcp), get_Df(qcp), dirty_quantifiers, dirty_qs, get_p(qcp), get_n(qcp), [X.v; intervals; get_G(qcp)])
        return test_zero_in(get_M(qcp), R_inner, get_positions_dict(qcp))
    end
end

function create_is_in_2(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where N
        quantifiers = [[(Exists, i) for i in 1:length(X)]; negation.(get_qvs(qcp)); [(Exists, get_p(qcp)-i) for i in (get_n(qcp)-1):-1:0]]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        qvs_relaxed = get_qvs_relaxed(qcp)
        quantifiers_relaxed = [[[(Exists, i) for i in 1:length(X)]; negation.(qvs_relaxed[j]); [(Exists, get_p(qcp)-i) for i in (get_n(qcp)-1):-1:0]] for j in 1:get_n(qcp)]
        dirty_qs = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        R_outer = QEapprox_o0_outer(get_f(qcp), get_Df(qcp), dirty_quantifiers, dirty_qs, get_p(qcp), get_n(qcp), [X.v; intervals; get_G(qcp)])
        return test_zero_not_in(get_M(qcp), R_outer, get_positions_dict(qcp)) 
    end
end

function create_is_out_1(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where N
        quantifiers = [[(Exists, i) for i in 1:length(X)]; get_qvs(qcp); [(Exists, get_p(qcp)-i) for i in (get_n(qcp)-1):-1:0]]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        qvs_relaxed = get_qvs_relaxed(qcp)
        quantifiers_relaxed = [[[(Exists, i) for i in 1:length(X)]; qvs_relaxed[j]; [(Exists, get_p(qcp)-i) for i in (get_n(qcp)-1):-1:0]] for j in 1:get_n(qcp)]
        dirty_qs = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        R_outer = QEapprox_o0_outer(get_f(qcp), get_Df(qcp), dirty_quantifiers, dirty_qs, get_p(qcp), get_n(qcp), [X.v; intervals; get_G(qcp)])
        return test_zero_not_in(get_M(qcp), R_outer, get_positions_dict(qcp))
    end
end

function create_is_out_2(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:Number}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where N
        quantifiers = [[(Forall, i) for i in 1:length(X)]; negation.(get_qvs(qcp)); [(Exists, get_p(qcp)-i) for i in (get_n(qcp)-1):-1:0]]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        qvs_relaxed = get_qvs_relaxed(qcp)
        quantifiers_relaxed = [[[(Forall, i) for i in 1:length(X)]; negation.(qvs_relaxed[j]); [(Exists, get_p(qcp)-i) for i in (get_n(qcp)-1):-1:0]] for j in 1:get_n(qcp)]
        dirty_qs = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        R_inner = QEapprox_o0_inner(get_f(qcp), get_Df(qcp), dirty_quantifiers, dirty_qs, get_p(qcp), get_n(qcp), [X.v; intervals; get_G(qcp)])
        return test_zero_in(get_M(qcp), R_inner, get_positions_dict(qcp))
    end
end

# function bounds(f, X, interval)
#     return [f[i]([X.v..., interval...]) for i in 1:length(f)]
# end

function check_is_in(X_0, P_in::Vector{Vector{IntervalArithmetic.Interval{T}}}, qcp, criterion) where T <: Number
    @assert criterion == 1 || criterion == 2

    indices_forall = [i for (q, i) in get_qvs(qcp) if q == Forall] .- length(X_0)
    indices_exists = [i for (q, i) in get_qvs(qcp) if q == Exists] .- length(X_0)

    indices = [1 for i in 1:length(P_in)]
    lengths = length.(P_in)
    is_in_union = false
    while !is_in_union && (isempty(indices_exists) || indices[indices_exists] <= lengths[indices_exists])
        is_in_intersection = true
        while is_in_intersection && (isempty(indices_forall) || indices[indices_forall] <= lengths[indices_forall])
            sub_interval = isempty(P_in) ? IntervalArithmetic.Interval{T}[] : [P_in[i][indices[i]] for i in 1:length(P_in)]
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

check_is_in_1(X_0, P_in, qcp) = check_is_in(X_0, P_in, qcp, 1)
check_is_in_2(X_0, P_in, qcp) = check_is_in(X_0, P_in, qcp, 2)

function check_is_out(X_0, P_out::Vector{Vector{IntervalArithmetic.Interval{T}}}, qcp, criterion) where T <: Number
    @assert criterion == 1 || criterion == 2

    indices_forall = [i for (q, i) in get_qvs(qcp) if q == Forall] .- length(X_0)
    indices_exists = [i for (q, i) in get_qvs(qcp) if q == Exists] .- length(X_0)

    indices = [1 for i in 1:length(P_out)]
    lengths = length.(P_out)
    is_out_intersection = true
    while is_out_intersection && (isempty(indices_exists) || indices[indices_exists] <= lengths[indices_exists])
        is_out_union = false
        while !is_out_union && (isempty(indices_forall) || indices[indices_forall] <= lengths[indices_forall])
            sub_interval = isempty(P_out) ? IntervalArithmetic.Interval{T}[] : [P_out[i][indices[i]] for i in 1:length(P_out)]
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

check_is_out_1(X_0, P_out, qcp) = check_is_out(X_0, P_out, qcp, 1)
check_is_out_2(X_0, P_out, qcp) = check_is_out(X_0, P_out, qcp, 2)

global z_in, z_out

const expansion_functions_1 = (expand_1, duplication_positions_1, complement_ranges_1)
const expansion_functions_2 = (expand_2, duplication_positions_2, complement_ranges_2)

function build_quantified_problem_in(parameters::ProblemParameters, domains, criterion)
    expand, duplication_positions, complement_ranges = criterion == 1 ? expansion_functions_1 : expansion_functions_2 
    
    M = get_M(parameters)
    userfunctions = get_f(parameters)
    qvs = get_qvs(parameters)
    qvs_relaxed = get_qvs_relaxed(parameters)
    sizes = get_sizes(parameters)
    n = get_n(parameters)
    p = get_p(parameters)
    
    G = get_G(domains)
    
    ranges_dict = build_ranges_dict(M, sizes)
    
    M_in = expand(M, ranges_dict, G)
    positions_in = duplication_positions(M, ranges_dict, G)
    Δn_in = length(positions_in)
    if userfunctions isa UserDefinedSymbolicFunctions
        variables = get_variables(userfunctions)
        f_num = get_f_num(userfunctions)
        @variables z_in[n + Δn_in]
        f_num_in = inflate(f_num, positions_in)
        for i in 1:(n + Δn_in)
            f_num_in[i] -= z_in[i]
        end
        f_fun_in, Df_fun_in = build_function_f_Df(f_num_in, [variables..., z_in...], n + Δn_in, p + n + Δn_in)
    elseif userfunctions isa UserDefinedFunctions
        f_inflated = inflate(get_f(userfunctions), positions_in)
        Df_inflated = inflate(get_Df(userfunctions), positions_in)
        f_fun_in = Function[]
        Df_fun_in = Function[]
        Dz = Matrix{Float64}(-LinearAlgebra.I, n + Δn_in, n + Δn_in)
        for i in 1:(n + Δn_in)
            f_fun_in_i = x -> f_inflated[i](x[1:p]) - x[p+i]
            push!(f_fun_in, f_fun_in_i)
            Df_fun_in_i = x -> [Df_inflated[i](x[1:p])..., Dz[i,:]...]
            push!(Df_fun_in, Df_fun_in_i)
        end
    else
        error("Unknown type of user-defined functions: $(typeof(userfunctions)).")
    end
    G_in = expand_intervals(G, complement_ranges(M, ranges_dict))
    qvs_relaxed_in = inflate(qvs_relaxed, positions_in)
    positions_dict_in = build_positions_dict(M_in)
    problem_in = ExpandedProblem(M_in, f_fun_in, Df_fun_in, G_in, positions_dict_in)
    
    return QuantifiedConstraintProblem(problem_in, qvs, qvs_relaxed_in, p + n + Δn_in, n + Δn_in)
end

function build_quantified_problem_out(parameters::ProblemParameters, domains, criterion)
    expand, duplication_positions, complement_ranges = criterion == 1 ? expansion_functions_1 : expansion_functions_2 
    
    M = get_M(parameters)
    userfunctions = get_f(parameters)
    qvs = get_qvs(parameters)
    qvs_relaxed = get_qvs_relaxed(parameters)
    sizes = get_sizes(parameters)
    n = get_n(parameters)
    p = get_p(parameters)
    
    G = get_G(domains)
    
    ranges_dict = build_ranges_dict(M, sizes)
    
    M_out = expand(M, ranges_dict, G)
    positions_out = duplication_positions(M, ranges_dict, G)
    Δn_out = length(positions_out)
    if userfunctions isa UserDefinedSymbolicFunctions
        variables = get_variables(userfunctions)
        f_num = get_f_num(userfunctions)
        @variables z_out[n + Δn_out]
        f_num_out = inflate(f_num, positions_out)
        for i in 1:(n + Δn_out)
            f_num_out[i] -= z_out[i]
        end
        f_fun_out, Df_fun_out = build_function_f_Df(f_num_out, [variables..., z_out...], n + Δn_out, p + n + Δn_out)
    elseif userfunctions isa UserDefinedFunctions
        f_inflated = inflate(get_f(userfunctions), positions_out)
        Df_inflated = inflate(get_Df(userfunctions), positions_out)
        f_fun_out = Function[]
        Df_fun_out = Function[]
        Dz = Matrix{Float64}(-LinearAlgebra.I, n + Δn_out, n + Δn_out)
        for i in 1:(n + Δn_out)
            f_fun_out_i = x -> f_inflated[i](x[1:p]) - x[p+i]
            push!(f_fun_out, f_fun_out_i)
            Df_fun_out_i = x -> [Df_inflated[i](x[1:p])..., Dz[i,:]...]
            push!(Df_fun_out, Df_fun_out_i)
        end
    else
        error("Unknown type of user-defined functions: $(typeof(userfunctions)).")
    end
    G_out = expand_intervals(G, complement_ranges(M, ranges_dict))
    qvs_relaxed_out = inflate(qvs_relaxed, positions_out)
    positions_dict_out = build_positions_dict(M_out)
    problem_out = ExpandedProblem(M_out, f_fun_out, Df_fun_out, G_out, positions_dict_out)
    
    return QuantifiedConstraintProblem(problem_out, qvs, qvs_relaxed_out, p + n + Δn_out, n + Δn_out)
end

Precision = Union{T, Vector{T}} where T<:Number

struct PavingConfiguration
    ϵ_x::Precision
    ϵ_p::Union{Precision, Nothing}
    allow_exists_or_forall_bisection::Bool
    allow_exists_and_forall_bisection::Bool
    function PavingConfiguration(ϵ_x::Precision)
        new(ϵ_x, nothing, false, false)
    end
    function PavingConfiguration(ϵ_x::Precision, ϵ_p::Precision, allow_exists_or_forall_bisection::Bool, allow_exists_and_forall_bisection::Bool)
        @assert nand(allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) "Refinement and subdivision are mutually exclusive."
        @assert ((allow_exists_and_forall_bisection || allow_exists_or_forall_bisection) && !isnothing(ϵ_p)) || (!allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection) "ϵ_p must be provided when bisection on parameter space is allowed."
        new(ϵ_x, ϵ_p, allow_exists_or_forall_bisection, allow_exists_and_forall_bisection)
    end
    function PavingConfiguration(ϵ_x::Precision, ϵ_p::Nothing, allow_exists_or_forall_bisection::Bool, allow_exists_and_forall_bisection::Bool)
        @assert !allow_exists_or_forall_bisection "ϵ_p must be provided when subdivision on parameter space is allowed."
        @assert !allow_exists_and_forall_bisection "ϵ_p must be provided when refinement on parameter space is allowed."
        new(ϵ_x, nothing, false, false)
    end
end

get_ϵ_x(configuration::PavingConfiguration) = configuration.ϵ_x
get_ϵ_p(configuration::PavingConfiguration) = configuration.ϵ_p
get_allow_exists_or_forall_bisection(configuration::PavingConfiguration) = configuration.allow_exists_or_forall_bisection
get_allow_exists_and_forall_bisection(configuration::PavingConfiguration) = configuration.allow_exists_and_forall_bisection
is_subdivided(configuration::PavingConfiguration) = configuration.allow_exists_or_forall_bisection
is_refined(configuration::PavingConfiguration) = configuration.allow_exists_and_forall_bisection

function Base.show(io::IO, configuration::PavingConfiguration)
    println(io, "ϵ_x: $(configuration.ϵ_x)")
    if !isnothing(configuration.ϵ_p)
        println(io, "ϵ_p: $(configuration.ϵ_p)")
    end
    println(io, if configuration.allow_exists_and_forall_bisection "Refined" else "Not refined" end)
    print(io, if configuration.allow_exists_or_forall_bisection "Normal bisection on P" else "No standard bisection on P" end)
end

function pave(X::IntervalArithmetic.IntervalBox{N, T}, parameters, domains, configuration, criterion_in, criterion_out)::Tuple{Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}} where {N, T<:Number}
    X_length = length(X)
    P_length = length(get_P(domains))
    @assert X_length + P_length == get_p(parameters) "Total number of variables, in X and p_in, must be equal to p = $(get_p(parameters))."
    @assert length(get_qvs(parameters)) == P_length "Number of quantified variables must be equal to the number of parameter."
    for qv in get_qvs(parameters)
        @assert X_length < index(qv) <= X_length + P_length "Quantified variables must be in the parameter space: indices between $(X_length+1) and $(X_length+P_length))."
    end
    for qvs in get_qvs_relaxed(parameters)
        @assert length(qvs) == length(get_P(domains)) "Number of quantified variables must be equal to the number of parameter boxes, $(P_length)."
        for qv in qvs
            @assert X_length < index(qv) <= X_length + P_length "Quantified variables must be in the parameter space: indices between $(X_length+1) and $(X_length+P_length)."
        end
    end
    check_is_in = criterion_in == 1 ? check_is_in_1 : check_is_in_2
    check_is_out = criterion_out == 1 ? check_is_out_1 : check_is_out_2
    
    n = get_n(parameters)
    p = get_p(parameters)
    
    qcp_in = build_quantified_problem_in(parameters, domains, criterion_in)
    qcp_out = build_quantified_problem_out(parameters, domains, criterion_out)
    
    P = get_P(domains)
    P_in = [[interval] for interval in P]
    P_out = deepcopy(P_in)

    ϵ_x = get_ϵ_x(configuration)
    ϵ_p = get_ϵ_p(configuration)
    allow_exists_or_forall_bisection = get_allow_exists_or_forall_bisection(configuration)
    allow_exists_and_forall_bisection = get_allow_exists_and_forall_bisection(configuration)

    if allow_exists_and_forall_bisection
        refine_in!(P_in, get_qvs(parameters), ϵ_p, p)
        refine_out!(P_out, get_qvs(parameters), ϵ_p, p)
    end

    P_in_0 = deepcopy(P_in)
    P_out_0 = deepcopy(P_out)

    inn = []
    out = []
    delta = []
    list = [(X, P_in, P_out)]
    while !isempty(list)
        X, P_in, P_out = pop!(list)
        if !allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection
            if check_is_in(X, P_in, qcp_in)
                push!(inn, X)
            elseif check_is_out(X, P_out, qcp_out)
                push!(out, X)
            elseif all(map(<, IntervalArithmetic.diam.(X), ϵ_x))
                push!(delta, X)
            else
                X_1, X_2, _ = bisect_precision(X, ϵ_x)
                push!(list, (X_1, deepcopy(P_in_0), deepcopy(P_out_0)))
                push!(list, (X_2, deepcopy(P_in_0), deepcopy(P_out_0)))
            end
        end
        if allow_exists_and_forall_bisection
            p_in_diams = IntervalArithmetic.diam.(first.(P_in))
            p_out_diams = IntervalArithmetic.diam.(first.(P_out))
            p_maxs = max.(p_in_diams, p_out_diams)
            X_diams = IntervalArithmetic.diam.(X)
            if check_is_in(X, P_in, qcp_in)
                push!(inn, X)
            elseif check_is_out(X, P_out, qcp_out)
                push!(out, X)
            elseif all(map(<, X_diams, ϵ_x)) && all(map(<, p_maxs, ϵ_p))
                push!(delta, X)
            elseif any(map(>=, X_diams, ϵ_x)) && all(map(<, p_maxs, ϵ_p))
                X_1, X_2, _ = bisect_precision(X, ϵ_x)
                push!(list, (X_1, deepcopy(P_in_0), deepcopy(P_out_0)))
                push!(list, (X_2, deepcopy(P_in_0), deepcopy(P_out_0)))
            elseif all(map(<, X_diams, ϵ_x)) && any(map(>=, p_maxs, ϵ_p))
                if any(map(<=, ϵ_p, p_in_diams))
                    bisect_largest_forall!(P_in, get_qvs(parameters), p, ϵ_p)
                    # bisect_eps_forall!(P_in, get_qvs(parameters), ϵ_p, p)
                end
                if any(map(<=, ϵ_p, p_out_diams))
                    bisect_largest_exists!(P_out, get_qvs(parameters), p, ϵ_p)
                    # bisect_eps_exists!(P_out, get_qvs(parameters), ϵ_p, p)
                end
                push!(list, (X, P_in, P_out))
            else
                if maximum(X_diams) < maximum(p_maxs)
                    if any(map(<=, ϵ_p, p_in_diams))
                        bisect_largest_forall!(P_in, get_qvs(parameters), p, ϵ_p)
                        # bisect_eps_forall!(P_in, get_qvs(parameters), ϵ_p, p)
                    end
                    if any(map(<=, ϵ_p, p_out_diams))
                        bisect_largest_exists!(P_out, get_qvs(parameters), p, ϵ_p)
                        # bisect_eps_exists!(P_out, get_qvs(parameters), ϵ_p, p)
                    end
                    push!(list, (X, P_in, P_out))
                else
                    X_1, X_2, _ = bisect_precision(X, ϵ_x)
                    push!(list, (X_1, deepcopy(P_in_0), deepcopy(P_out_0)))
                    push!(list, (X_2, deepcopy(P_in_0), deepcopy(P_out_0)))
                end
            end
        end
        if allow_exists_or_forall_bisection
            indices_forall = [i for (q, i) in qvs if q == Forall] .- length(X)
            indices_exists = [i for (q, i) in qvs if q == Exists] .- length(X)
            p_in_diams = IntervalArithmetic.diam.(first.(P_in))
            p_in_diams[indices_forall] .= -1.0
            p_out_diams = IntervalArithmetic.diam.(first.(P_out))
            p_out_diams[indices_exists] .= -1.0
            p_maxs = max.(p_in_diams, p_out_diams)
            X_diams = IntervalArithmetic.diam.(X)
            if check_is_in(X, P_in, qcp_in)
                push!(inn, X)
            elseif check_is_out(X, P_out, qcp_out)
                push!(out, X)
            elseif all(map(<, X_diams, ϵ_x)) && all(map(<, p_maxs, ϵ_p))
                push!(delta, X)
            elseif any(map(>=, X_diams, ϵ_x)) && all(map(<, p_maxs, ϵ_p))
                X_1, X_2, _ = bisect_precision(X, ϵ_x)
                push!(list, (X_1, deepcopy(P_in_0), deepcopy(P_out_0)))
                push!(list, (X_2, deepcopy(P_in_0), deepcopy(P_out_0)))
            elseif all(map(<, X_diams, ϵ_x)) && any(map(>=, p_maxs, ϵ_p))
                if any(map(<=, ϵ_p, p_in_diams))
                    bisect_largest_forall!(P_in, get_qvs(qcp_in), p, ϵ_p)
                end
                if any(map(<=, ϵ_p, p_out_diams))
                    bisect_largest_exists!(P_out, get_qvs(qcp_out), p, ϵ_p)
                end
                push!(list, (X, P_in, P_out))
            else
                if maximum(X_diams) < maximum(p_maxs)
                    if any(map(<=, ϵ_p, p_in_diams))
                        bisect_largest_exists!(P_in, get_qvs(qcp_in), p, ϵ_p)
                    end
                    if any(map(<=, ϵ_p, p_out_diams))
                        bisect_largest_forall!(P_out, get_qvs(qcp_out), p, ϵ_p)
                    end
                    push!(list, (X, P_in, P_out))
                else
                    X_1, X_2, _ = bisect_precision(X, ϵ_x)
                    push!(list, (X_1, deepcopy(P_in_0), deepcopy(P_out_0)))
                    push!(list, (X_2, deepcopy(P_in_0), deepcopy(P_out_0)))
                end
            end
        end
    end
    return inn, out, delta
end

pave_11(X, parameters, domains, configuration) = pave(X, parameters, domains, configuration, 1, 1)
pave_12(X, parameters, domains, configuration) = pave(X, parameters, domains, configuration, 1, 2)
pave_21(X, parameters, domains, configuration) = pave(X, parameters, domains, configuration, 2, 1)
pave_22(X, parameters, domains, configuration) = pave(X, parameters, domains, configuration, 2, 2)

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

function backward_bound(interval, sign)
    if sign > 0
        return interval.lo
    else
        return interval.hi
    end
end

function forward_bound(interval, sign)
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

backward_slice(box, optimization_directions, dim) = slice(box, optimization_directions, dim, backward_bound)
forward_slice(box, optimization_directions, dim) = slice(box, optimization_directions, dim, forward_bound)

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