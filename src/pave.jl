const global minus_inf = -100000
const global plus_inf = 100000
const global strict_epsilon = 0.0001
global NB_EVALUATIONS::Unsigned

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

function bisect_eps_quantifier!(intervals, qvs, eps, shift, quantifier)
    @assert sum(length.(intervals); init=0) ==  length(intervals) "Each interval should be a single interval."
    pos_quantifier = [i for (q, i) in qvs if q == quantifier] .- shift

    for i in pos_quantifier
        intervals[i] = bisect_eps(intervals[i][1], eps[i])
    end
end

bisect_eps_exists!(intervals, qvs, eps, shift) = bisect_eps_quantifier!(intervals, qvs, eps, shift, Exists)
bisect_eps_forall!(intervals, qvs, eps, shift) = bisect_eps_quantifier!(intervals, qvs, eps, shift, Forall)

function pointify_quantifier!(intervals, qvs, shift, quantifier)
    pos_quantifier = [i for (q, i) in qvs if q == quantifier] .- shift

    for i in pos_quantifier
        intervals[i] = interval.(mid.(intervals[i]))
    end
end

pointify_exists!(intervals, qvs, shift) = pointify_quantifier!(intervals, qvs, shift, Exists)
pointify_forall!(intervals, qvs, shift) = pointify_quantifier!(intervals, qvs, shift, Forall)

function refine_in!(P_in, qvs, eps, shift)
    bisect_eps_exists!(P_in, qvs, eps, shift)
    pointify_exists!(P_in, qvs, shift)
end

function refine_out!(P_out, qvs, eps, shift)
    bisect_eps_forall!(P_out, qvs, eps, shift)
    pointify_forall!(P_out, qvs, shift)
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

function complement(x::IntervalArithmetic.Interval{T}) where T <: AbstractFloat
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
function is_bounded(interval::IntervalArithmetic.Interval{T}) where T <: AbstractFloat
    return !isempty(interval) && interval.lo != minus_inf && interval.hi != plus_inf
end

function is_bounded(intervals::Vector{IntervalArithmetic.Interval{T}}) where T <: AbstractFloat
    return all(is_bounded, intervals)
end

function atom_conjunction(atom::Atom, n::Int)
    if n == 1
        return atom
    end
    conjuncted_predicates = [Atom("$(atom.value)_$(i)") for i in 1:n]
    return ∧(conjuncted_predicates...)
end

function atom_complement_disjunction(atom::Atom, intervals::Vector{IntervalArithmetic.Interval{T}}) where T <: AbstractFloat
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

function expand_1(M, indices_dict, intervals::Vector{IntervalArithmetic.Interval{T}}, atoms) where T <: AbstractFloat
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

function expand_2(M, indices_dict, intervals::Vector{IntervalArithmetic.Interval{T}}, atoms, previous_token=nothing) where T <: AbstractFloat
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

function expand_intervals(intervals::Vector{IntervalArithmetic.Interval{T}}, ranges::Vector{UnitRange{Int}}) where T <: AbstractFloat
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

function inflate(vector::AbstractVector, indices)
    inflated = []
    for (i, item) in enumerate(vector)
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

function create_is_in_1(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:AbstractFloat}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where N
        qvs_relaxed = get_qvs_relaxed(qcp)
        quantifiers_relaxed = [[[(Forall, i) for i in 1:length(X)]; qvs_relaxed[j]; [(Exists, length(X)+length(qvs_relaxed[j])+i) for i in 1:length(get_G(qcp))]] for j in eachindex(qvs_relaxed)]
        dirty_quantifiers_relaxed = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        R_inner = QEapprox_o0_inner(get_f(qcp), get_Df(qcp), dirty_quantifiers_relaxed, length(get_G(qcp)), [X.v; intervals; get_G(qcp)])
        global NB_EVALUATIONS
        NB_EVALUATIONS += 1
        return test_zero_in(get_M(qcp), R_inner, get_positions_dict(qcp))
    end
end

function create_is_in_2(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:AbstractFloat}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where N
        quantifiers = [[(Exists, i) for i in 1:length(X)]; negation.(get_qvs(qcp)); [(Exists, length(X)+length(get_qvs(qcp))+i) for i in 1:length(get_G(qcp))]]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        R_outer = QEapprox_o0_outer(get_f(qcp), get_Df(qcp), dirty_quantifiers, length(get_G(qcp)), [X.v; intervals; get_G(qcp)])
        global NB_EVALUATIONS
        NB_EVALUATIONS += 1
        return test_zero_not_in(get_M(qcp), R_outer, get_positions_dict(qcp)) 
    end
end

function create_is_out_1(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:AbstractFloat}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where N
        quantifiers = [[(Exists, i) for i in 1:length(X)]; get_qvs(qcp); [(Exists, length(X)+length(get_qvs(qcp))+i) for i in 1:length(get_G(qcp))]]
        dirty_quantifiers = quantifiedvariables2dirtyvariables(quantifiers)
        R_outer = QEapprox_o0_outer(get_f(qcp), get_Df(qcp), dirty_quantifiers, length(get_G(qcp)), [X.v; intervals; get_G(qcp)])
        global NB_EVALUATIONS
        NB_EVALUATIONS += 1
        return test_zero_not_in(get_M(qcp), R_outer, get_positions_dict(qcp))
    end
end

function create_is_out_2(qcp::QuantifiedConstraintProblem, intervals::AbstractVector{IntervalArithmetic.Interval{T}})::Function where {T<:AbstractFloat}
    return function(X::IntervalArithmetic.IntervalBox{N, T}) where N
        qvs_relaxed = get_qvs_relaxed(qcp)
        quantifiers_relaxed = [[[(Forall, i) for i in 1:length(X)]; negation.(qvs_relaxed[j]); [(Exists, length(X)+length(qvs_relaxed[j])+i) for i in 1:length(get_G(qcp))]] for j in eachindex(qvs_relaxed)]
        dirty_quantifiers_relaxed = quantifiedvariables2dirtyvariables.(quantifiers_relaxed)
        R_inner = QEapprox_o0_inner(get_f(qcp), get_Df(qcp), dirty_quantifiers_relaxed, length(get_G(qcp)), [X.v; intervals; get_G(qcp)])
        global NB_EVALUATIONS
        NB_EVALUATIONS += 1
        return test_zero_in(get_M(qcp), R_inner, get_positions_dict(qcp))
    end
end

# function bounds(f, X, interval)
#     return [f[i]([X.v..., interval...]) for i in 1:length(f)]
# end

function check_is_in(X_0, P_in::Vector{Vector{IntervalArithmetic.Interval{T}}}, qcp, criterion) where T <: AbstractFloat
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

function check_is_out(X_0, P_out::Vector{Vector{IntervalArithmetic.Interval{T}}}, qcp, criterion) where T <: AbstractFloat
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
    p = get_p(parameters)
    
    G = get_G(domains)
    
    ranges_dict = build_ranges_dict(M, sizes)
    
    M_in = expand(M, ranges_dict, G)
    n_in = length(atoms(M_in))
    positions_in = duplication_positions(M, ranges_dict, G)
    # Δn_in = length(positions_in)
    if userfunctions isa UserDefinedSymbolicFunctions
        variables = get_variables(userfunctions)
        f_num = get_f_num(userfunctions)
        @variables z_in[1:n_in]
        f_num_in = inflate(f_num, positions_in)
        for i in 1:(n_in)
            f_num_in[i] -= z_in[i]
        end
        f_fun_in, Df_fun_in = build_function_f_Df_2(f_num_in, [variables; z_in], n_in, p + n_in)
        Df_fun_in = MatJacobianMatrix(Df_fun_in)
    elseif userfunctions isa UserDefinedFunctions
        f_inflated = inflate(get_f(userfunctions), positions_in)
        Df_inflated = inflate(get_Df(userfunctions), positions_in)
        f_fun_in = Function[]
        # Df_fun_in = Matrix{Function}(undef, n_in, p + n_in)
        Dz = Matrix{Float64}(-LinearAlgebra.I, n_in, n_in)
        for i in 1:(n_in)
            f_fun_in_i = (x...) -> f_inflated[i](x[1:p]...) - x[p+i]
            push!(f_fun_in, f_fun_in_i)
        end
        Df_fun_in = append_constant(Df_inflated, Dz)
    else
        error("Unknown type of user-defined functions: $(typeof(userfunctions)).")
    end
    G_in = expand_intervals(G, complement_ranges(M, ranges_dict))
    qvs_relaxed_in = inflate(qvs_relaxed, positions_in)
    positions_dict_in = build_positions_dict(M_in)
    problem_in = ExpandedProblem(M_in, f_fun_in, Df_fun_in, G_in, positions_dict_in)
    
    return QuantifiedConstraintProblem(problem_in, qvs, qvs_relaxed_in)
end

function build_quantified_problem_out(parameters::ProblemParameters, domains, criterion)
    expand, duplication_positions, complement_ranges = criterion == 1 ? expansion_functions_1 : expansion_functions_2 
    
    M = get_M(parameters)
    userfunctions = get_f(parameters)
    qvs = get_qvs(parameters)
    qvs_relaxed = get_qvs_relaxed(parameters)
    sizes = get_sizes(parameters)
    p = get_p(parameters)
    
    G = get_G(domains)
    
    ranges_dict = build_ranges_dict(M, sizes)
    
    M_out = expand(M, ranges_dict, G)
    n_out = length(atoms(M_out))
    positions_out = duplication_positions(M, ranges_dict, G)
    # Δn_out = length(positions_out)
    if userfunctions isa UserDefinedSymbolicFunctions
        variables = get_variables(userfunctions)
        f_num = get_f_num(userfunctions)
        @variables z_out[1:n_out]
        f_num_out = inflate(f_num, positions_out)
        for i in 1:(n_out)
            f_num_out[i] -= z_out[i]
        end
        f_fun_out, Df_fun_out = build_function_f_Df_2(f_num_out, [variables; z_out], n_out, p + n_out)
        Df_fun_out = MatJacobianMatrix(Df_fun_out)
    elseif userfunctions isa UserDefinedFunctions
        f_inflated = inflate(get_f(userfunctions), positions_out)
        Df_inflated = inflate(get_Df(userfunctions), positions_out)
        f_fun_out = Function[]
        # Df_fun_out = Matrix{Function}(undef, n_out, p + n_out)
        Dz = Matrix{Float64}(-LinearAlgebra.I, n_out, n_out)
        for i in 1:(n_out)
            f_fun_out_i = (x...) -> f_inflated[i](x[1:p]...) - x[p+i]
            push!(f_fun_out, f_fun_out_i)
        end
        Df_fun_out = append_constant(Df_inflated, Dz)
    else
        error("Unknown type of user-defined functions: $(typeof(userfunctions)).")
    end
    G_out = expand_intervals(G, complement_ranges(M, ranges_dict))
    qvs_relaxed_out = inflate(qvs_relaxed, positions_out)
    positions_dict_out = build_positions_dict(M_out)
    problem_out = ExpandedProblem(M_out, f_fun_out, Df_fun_out, G_out, positions_dict_out)
    
    return QuantifiedConstraintProblem(problem_out, qvs, qvs_relaxed_out)
end

struct PavingConfiguration{T <: AbstractFloat}
    ϵ_x::AbstractVector{T}
    ϵ_p::AbstractVector{T}
    allow_exists_or_forall_bisection::Bool
    allow_exists_and_forall_bisection::Bool
    function PavingConfiguration(ϵ_x)
        PavingConfiguration(ϵ_x, [], false, false)
    end
    function PavingConfiguration(ϵ_x, ϵ_p, allow_exists_or_forall_bisection, allow_exists_and_forall_bisection)
        type_x = ϵ_x isa AbstractVector ? eltype(ϵ_x) : typeof(ϵ_x)
        @assert type_x <: AbstractFloat "ϵ_x must be a floating-point number or a vector of floating-point numbers."
        PavingConfiguration(isempty(ϵ_x) ? type_x[] : float(ϵ_x), isempty(ϵ_p) ? type_x[] : float(ϵ_p), allow_exists_or_forall_bisection, allow_exists_and_forall_bisection)
    end
    function PavingConfiguration(ϵ_x::T, ϵ_p, allow_exists_or_forall_bisection, allow_exists_and_forall_bisection) where T <: AbstractFloat
        PavingConfiguration([ϵ_x], ϵ_p, allow_exists_or_forall_bisection, allow_exists_and_forall_bisection)
    end
    function PavingConfiguration(ϵ_x::AbstractVector{T}, ϵ_p::T, allow_exists_or_forall_bisection::Bool, allow_exists_and_forall_bisection::Bool) where T <: AbstractFloat
        PavingConfiguration(ϵ_x, [ϵ_p], allow_exists_or_forall_bisection, allow_exists_and_forall_bisection)
    end
    function PavingConfiguration(ϵ_x::AbstractVector{T}, ϵ_p::AbstractVector{T}, allow_exists_or_forall_bisection::Bool, allow_exists_and_forall_bisection::Bool) where T <: AbstractFloat
        @assert !isempty(ϵ_x) "ϵ_x must be provided."
        @assert nand(allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) "Refinement and subdivision are mutually exclusive."
        @assert ((allow_exists_and_forall_bisection || allow_exists_or_forall_bisection) && !isempty(ϵ_p)) || (!allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection) "ϵ_p must be provided when bisection on parameter space is allowed."
        new{T}(ϵ_x, ϵ_p, allow_exists_or_forall_bisection, allow_exists_and_forall_bisection)
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

function pave(X::IntervalArithmetic.IntervalBox{N, T}, parameters, domains, configuration, criterion_in, criterion_out, verbose=true)::Tuple{Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}, Vector{IntervalArithmetic.IntervalBox{N, T}}} where {N, T<:AbstractFloat}
    X_length = length(X)
    P_length = length(get_P(domains))
    @assert X_length + P_length == get_p(parameters) "Total number of variables, in X and P, must be equal to p = $(get_p(parameters))."
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
    
    # n = get_n(parameters)
    p = get_p(parameters)
    
    qcp_in = build_quantified_problem_in(parameters, domains, criterion_in)
    qcp_out = build_quantified_problem_out(parameters, domains, criterion_out)
    
    P = get_P(domains)
    P_in = [[interval] for interval in P]
    P_out = deepcopy(P_in)

    ϵ_x = get_ϵ_x(configuration)
    @assert length(ϵ_x) == X_length "Length of ϵ_x must be equal to the number of variables in X: $(X_length)."
    ϵ_p = get_ϵ_p(configuration)
    @assert isempty(ϵ_p) || length(ϵ_p) == P_length "Length of ϵ_p must be equal to the number of parameters: $(P_length)."
    allow_exists_or_forall_bisection = get_allow_exists_or_forall_bisection(configuration)
    allow_exists_and_forall_bisection = get_allow_exists_and_forall_bisection(configuration)

    if allow_exists_and_forall_bisection
        refine_in!(P_in, get_qvs(parameters), ϵ_p, X_length)
        refine_out!(P_out, get_qvs(parameters), ϵ_p, X_length)
    end

    P_in_0 = deepcopy(P_in)
    P_out_0 = deepcopy(P_out)

    println("Paving of the following domain")
    display(X)
    println("Target domain")
    display(get_G(domains))
    inn = []
    out = []
    delta = []
    list = [(X, P_in, P_out)]
    nb_iterations = 0
    global NB_EVALUATIONS
    NB_EVALUATIONS = 0
    while !isempty(list)
        nb_iterations += 1
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
    if verbose
        println("Paving completed.")
        println("\tNumber of iterations: $nb_iterations")
        println("\tNumber of evaluations: $NB_EVALUATIONS")
    end
    return inn, out, delta
end

pave(X::AbstractVector{IntervalArithmetic.Interval{T}}, parameters, domains, configuration, verbose=true) where T <: AbstractFloat = pave(IntervalBox(X), parameters, domains, configuration, verbose)
pave(X::IntervalArithmetic.Interval{T}, parameters, domains, configuration, verbose=true) where T <: AbstractFloat = pave([X], parameters, domains, configuration, verbose)

pave_11(X, parameters, domains, configuration, verbose=true) = pave(X, parameters, domains, configuration, 1, 1, verbose)
pave_12(X, parameters, domains, configuration, verbose=true) = pave(X, parameters, domains, configuration, 1, 2, verbose)
pave_21(X, parameters, domains, configuration, verbose=true) = pave(X, parameters, domains, configuration, 2, 1, verbose)
pave_22(X, parameters, domains, configuration, verbose=true) = pave(X, parameters, domains, configuration, 2, 2, verbose)

pave(X, parameters, domains, configuration, verbose=true) = pave_21(X, parameters, domains, configuration, verbose)

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
function pave_monotonous_mid(X::IntervalArithmetic.IntervalBox{N, T}, optimization_directions, parameters, domains, configuration, criterion_in, criterion_out, verbose=true) where {N, T<:AbstractFloat}
    @assert length(optimization_directions) == length(X) "Length of optimization_directions must be equal to the number of variables in X."
    X_length = length(X)
    P_length = length(get_P(domains))
    @assert X_length + P_length == get_p(parameters) "Total number of variables, in X and P, must be equal to p = $(get_p(parameters))."
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
    
    p = get_p(parameters)
    
    qcp_in = build_quantified_problem_in(parameters, domains, criterion_in)
    qcp_out = build_quantified_problem_out(parameters, domains, criterion_out)
    
    P = get_P(domains)
    P_in = [[interval] for interval in P]
    P_out = deepcopy(P_in)

    ϵ_x = get_ϵ_x(configuration)
    @assert length(ϵ_x) == X_length "Length of ϵ_x must be equal to the number of variables in X: $(X_length)."
    ϵ_p = get_ϵ_p(configuration)
    @assert isempty(ϵ_p) || length(ϵ_p) == P_length "Length of ϵ_p must be equal to the number of parameters: $(P_length)."
    allow_exists_or_forall_bisection = get_allow_exists_or_forall_bisection(configuration)
    allow_exists_and_forall_bisection = get_allow_exists_and_forall_bisection(configuration)

    # P_in_0 = deepcopy(P_in)
    # P_out_0 = deepcopy(P_out)
    P_in_diams = IntervalArithmetic.diam.(first.(P_in))
    P_out_diams = IntervalArithmetic.diam.(first.(P_out))
    if any(map(<=, ϵ_p, P_in_diams))
        # bisect_largest_forall!(p_in, qcp.qvs, qcp.p, qcp.n)
        bisect_eps_forall!(P_in, get_qvs(qcp_in), ϵ_p, X_length)
    end
    if any(map(<=, ϵ_p, P_out_diams))
        # bisect_largest_exists!(p_out, qcp.qvs, qcp.p, qcp.n)
        bisect_eps_exists!(P_out, get_qvs(qcp_out), ϵ_p, X_length)
    end

    inn = []
    out = []
    delta = []
    stack = [(X, P_in, P_out)]
    global NB_EVALUATIONS
    NB_EVALUATIONS = 0
    while !isempty(stack)
        X, P_in, P_out = pop!(stack)
        if !allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection
            error("TO DO")
        end
        if allow_exists_and_forall_bisection
            if all(map(<, IntervalArithmetic.diam.(X), ϵ_x))
                push!(delta, X)
                continue
            end
            dim_slice, X_slice = bisection_slice(X, ϵ_x)
            if check_is_in(X_slice, P_in, qcp_in)
                bisect_in_optimization_direction = create_bisect_in_optimization_direction(optimization_directions[dim_slice])
                X_inn, X_stack = bisect_in_optimization_direction(X, dim_slice)
                push!(inn, X_inn)
                push!(stack, (X_stack, P_in, P_out))
            elseif check_is_out(X_slice, P_out, qcp_out)
                bisect_in_optimization_direction = create_bisect_in_optimization_direction(optimization_directions[dim_slice])
                X_stack, X_out = bisect_in_optimization_direction(X, dim_slice)
                push!(out, X_out)
                push!(stack, (X_stack, P_in, P_out))
            else
                X_1, X_2, _ = bisect_precision(X, ϵ_x)
                push!(stack, (X_1, P_in, P_out))
                push!(stack, (X_2, P_in, P_out))
            end
        end
        if allow_exists_or_forall_bisection
            error("TO DO")
        end
    end
    if verbose
        println("Paving completed.")
        println("\tNumber of evaluations: $NB_EVALUATIONS")
    end
    return inn, out, delta
end

pave_monotonous_mid_12(X, optimization_directions, parameters, domains, configuration, verbose=true) = pave_monotonous_mid(X, optimization_directions, parameters, domains, configuration, 1, 2, verbose)

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
        if !isempty(unchanged_face) && i == dim(unchanged_face) && bound(directed_slicing_function, optimization_directions[i]) == bound(unchanged_face)
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

struct Face
    dim::Int
    bound::Bound
end

dim(face::Face) = face.dim
bound(face::Face) = face.bound

emptyface() = Face(-1, Upper)
Base.isempty(face::Face) = face.dim == -1

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
function pave_monotonous_sides_optimized(X::IntervalArithmetic.IntervalBox{N, T}, optimization_directions, parameters, domains, configuration, criterion_in, criterion_out, verbose=true) where {N, T<:AbstractFloat}
    @assert length(optimization_directions) == length(X) "Length of optimization_directions must be equal to the number of variables in X."
    X_length = length(X)
    P_length = length(get_P(domains))
    @assert X_length + P_length == get_p(parameters) "Total number of variables, in X and P, must be equal to p = $(get_p(parameters))."
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
    
    p = get_p(parameters)
    
    qcp_in = build_quantified_problem_in(parameters, domains, criterion_in)
    qcp_out = build_quantified_problem_out(parameters, domains, criterion_out)
    
    P = get_P(domains)
    P_in = [[interval] for interval in P]
    P_out = deepcopy(P_in)
    
    ϵ_x = get_ϵ_x(configuration)
    @assert length(ϵ_x) == X_length "Length of ϵ_x must be equal to the number of variables in X: $(X_length)."
    ϵ_p = get_ϵ_p(configuration)
    @assert isempty(ϵ_p) || length(ϵ_p) == P_length "Length of ϵ_p must be equal to the number of parameters: $(P_length)."
    allow_exists_or_forall_bisection = get_allow_exists_or_forall_bisection(configuration)
    allow_exists_and_forall_bisection = get_allow_exists_and_forall_bisection(configuration)
    
    # P_in_0 = deepcopy(P_in)
    # P_out_0 = deepcopy(P_out)
    P_in_diams = IntervalArithmetic.diam.(first.(P_in))
    P_out_diams = IntervalArithmetic.diam.(first.(P_out))
    if any(map(<=, ϵ_p, P_in_diams))
        # bisect_largest_forall!(p_in, qcp.qvs, qcp.p, qcp.n)
        bisect_eps_forall!(P_in, get_qvs(qcp_in), ϵ_p, X_length)
    end
    if any(map(<=, ϵ_p, P_out_diams))
        # bisect_largest_exists!(p_out, qcp.qvs, qcp.p, qcp.n)
        bisect_eps_exists!(P_out, get_qvs(qcp_out), ϵ_p, X_length)
    end
    inn = []
    out = []
    delta = []
    queue = [(X, emptyface(), P_in, P_out)]
    global NB_EVALUATIONS
    NB_EVALUATIONS = 0
    while !isempty(queue)
        X, unchanged_face, P_in, P_out = popfirst!(queue)
        if !allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection
            error("TO DO")
        end
        if allow_exists_and_forall_bisection
            backward_X_slices = backward_slices(X, optimization_directions, unchanged_face)
            forward_X_slices = forward_slices(X, optimization_directions, unchanged_face)
            if any(check_is_in(X_slice, P_in, qcp_in) for X_slice in forward_X_slices)
                push!(inn, X)
            elseif any(check_is_out(X_slice, P_out, qcp_out) for X_slice in backward_X_slices)
                push!(out, X)
            elseif all(map(<, IntervalArithmetic.diam.(X), ϵ_x))
                push!(delta, X)
            else
                X_1, X_2, dim = bisect_precision(X, ϵ_x)
                unchanged_face_1 = Face(dim, Lower)
                unchanged_face_2 = Face(dim, Upper)
                push!(queue, (X_1, unchanged_face_1, P_in, P_out))
                push!(queue, (X_2, unchanged_face_2, P_in, P_out))
            end
        end
        if allow_exists_or_forall_bisection
            error("TO DO")
        end
    end
    if verbose
        println("Paving completed.")
        println("\tNumber of evaluations: $NB_EVALUATIONS")
    end
    return inn, out, delta
end

pave_monotonous_sides_optimized_12(X, optimization_directions, parameters, domains, configuration, verbose=true) = pave_monotonous_sides_optimized(X, optimization_directions, parameters, domains, configuration, 1, 2, verbose)
