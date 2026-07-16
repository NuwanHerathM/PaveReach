#!/usr/local/bin/julia

using Match
using ReusePatterns
using SoleLogics

@enum Quantifier begin
    Forall
    Exists
end

QuantifiedVariable = Tuple{Quantifier, Int}

function quantifier(qv::QuantifiedVariable)
    return qv[1]
end

function index(qv::QuantifiedVariable)
    return qv[2]
end

function negation(quantifier::Quantifier)
    @match quantifier begin
        $Forall => return Exists
        $Exists => return Forall
        _ => error("Unknown quantifier: $quantifier")
    end
end

function negation(qv::QuantifiedVariable)
    quantifier, index = qv
    return (negation(quantifier), index)
end

function quantifiedvariables2dirtyvariables(qvs::Vector{QuantifiedVariable})
    dirty_variables = Any[]
    for qv in qvs
        quantifier, index = qv
        @match quantifier begin
            $Forall => push!(dirty_variables, "forall", index)
            $Exists => push!(dirty_variables, "exists", index)
            _ => error("Unknown quantifier: $quantifier")
        end
    end
    return dirty_variables
end

#----------------------------------------------------------------------------------------

function build_ranges(sizes::Vector{Int})
    end_indices::Vector{Int} = cumsum(sizes)
    start_indices::Vector{Int} = end_indices .- sizes .+ 1
    return range.(start_indices, end_indices)
end

function build_ranges_dict(M::SyntaxTree, sizes::Vector{Int})
    keys = map(atom -> atom.value, atoms(M))
    ranges = build_ranges(sizes)
    return Dict(zip(keys, ranges))
end

function build_positions_dict(M::SyntaxTree)
    atom_from_position = Dict(enumerate(atoms(M)))
    position_from_atom = Dict(value => key for (key, value) in atom_from_position)
    return position_from_atom
end

SymbolicFunction = Num

struct UserDefinedFunctions
    f::Vector{Function}
    Df::Vector{Function}
end

get_f(userfunctions::UserDefinedFunctions) = userfunctions.f
get_Df(userfunctions::UserDefinedFunctions) = userfunctions.Df

struct UserDefinedSymbolicFunctions
    variables::AbstractVector{Num}
    f_num::Vector{Num}
end

get_variables(userfunctions::UserDefinedSymbolicFunctions) = userfunctions.variables
get_f_num(userfunctions::UserDefinedSymbolicFunctions) = userfunctions.f_num

struct ProblemParameters
    M::SyntaxTree
    f::Union{UserDefinedFunctions, UserDefinedSymbolicFunctions}
    sizes::Vector{Int}
    qvs::Vector{QuantifiedVariable}
    qvs_relaxed::Vector{Vector{QuantifiedVariable}}
    n::Int
    p::Int
    function ProblemParameters(variables::AbstractVector{Num}, f_num::Vector{Num}, qvs::Vector{QuantifiedVariable}, qvs_relaxed::Vector{Vector{QuantifiedVariable}}, p::Int)
        n = length(f_num)
        M = parseformula("P")
        sizes = [n]
        ProblemParameters(M, variables, f_num, sizes, qvs, qvs_relaxed, p)
    end
    function ProblemParameters(M::SyntaxTree, variables::AbstractVector{Num}, f_num::Vector{Num}, sizes::Vector{Int}, qvs::Vector{QuantifiedVariable}, qvs_relaxed::Vector{Vector{QuantifiedVariable}}, p::Int)
        n = length(f_num)
        @assert length(qvs_relaxed) == n "Number of relaxed quantifier variable lists should be equal to number of functions."
        @assert length(sizes) == length(atoms(M)) "Number of dimensions should be equal to number of predicates."
        userfunctions = UserDefinedSymbolicFunctions(variables, f_num)
        new(M, userfunctions, sizes, qvs, qvs_relaxed, n, p)
    end
    function ProblemParameters(M::SyntaxTree, f::Vector{Function}, Df::Vector{Function}, sizes::Vector{Int}, qvs::Vector{QuantifiedVariable}, qvs_relaxed::Vector{Vector{QuantifiedVariable}}, p::Int)
        n = length(f)
        @assert length(Df) == n "Number of functions and gradients must be equal."
        @assert length(qvs_relaxed) == n "Number of relaxed quantifier variable lists should be equal to number of functions."
        @assert length(sizes) == length(atoms(M)) "Number of dimensions should be equal to number of predicates."
        userfunctions = UserDefinedFunctions(f, Df)
        new(M, userfunctions, sizes, qvs, qvs_relaxed, n, p)
    end
end

get_M(parameters::ProblemParameters) = parameters.M
get_f(parameters::ProblemParameters) = parameters.f
get_qvs(parameters::ProblemParameters) = parameters.qvs
get_qvs_relaxed(parameters::ProblemParameters) = parameters.qvs_relaxed
get_sizes(parameters::ProblemParameters) = parameters.sizes
get_n(parameters::ProblemParameters) = parameters.n
get_p(parameters::ProblemParameters) = parameters.p

struct ProblemDomains{T}
    P::Vector{IntervalArithmetic.Interval{T}}
    G::Vector{IntervalArithmetic.Interval{T}}
end

get_P(domains::ProblemDomains) = domains.P
get_G(domains::ProblemDomains) = domains.G

struct ExpandedProblem{T}
    M::SyntaxTree
    f::Vector{Function}
    Df::Vector{Function}
    G::Vector{IntervalArithmetic.Interval{T}}
    positions_dict::Dict{Atom{String}, Int}
    function ExpandedProblem(M::SyntaxTree, f::Vector{Function}, Df::Vector{Function}, G::Vector{IntervalArithmetic.Interval{T}}, positions_dict::Dict{Atom{String}, Int}) where T <: Number
        @assert allunique(atoms(M)) "All predicates in the formula must be unique."
        @assert length(f) == length(Df) "Number of functions and gradients must be equal."
        @assert length(f) == length(G) "Number of functions and domains must be equal."
        @assert length(f) == length(positions_dict) "Number of functions and positions must be equal."
        @assert length(positions_dict) == length(atoms(M)) "Number of positions and predicates must be equal."
        new{T}(M, f, Df, G, positions_dict)
    end
end

get_M(problem::ExpandedProblem) = problem.M
get_f(problem::ExpandedProblem) = problem.f
get_Df(problem::ExpandedProblem) = problem.Df
get_G(problem::ExpandedProblem) = problem.G
get_positions_dict(problem::ExpandedProblem) = problem.positions_dict

struct QuantifiedConstraintProblem
    problem::ExpandedProblem
    qvs::Vector{QuantifiedVariable}
    qvs_relaxed::Vector{Vector{QuantifiedVariable}}
    p::Int
    n::Int
    function QuantifiedConstraintProblem(problem::ExpandedProblem, qvs, qvs_relaxed, p::Int, n::Int)
        @assert length(problem.f) == n "Number of functions in the problem should be equal to n = $n."
        @assert length(qvs_relaxed) == n "Number of relaxed quantifier variable lists should be equal to n = $n."
        new(problem, qvs, qvs_relaxed, p, n)
    end
end

@forward((QuantifiedConstraintProblem, :problem), ExpandedProblem)
get_qvs(qcp::QuantifiedConstraintProblem) = qcp.qvs
get_qvs_relaxed(qcp::QuantifiedConstraintProblem) = qcp.qvs_relaxed
get_n(qcp::QuantifiedConstraintProblem) = qcp.n
get_p(qcp::QuantifiedConstraintProblem) = qcp.p

function quantifier(qcp::QuantifiedConstraintProblem, dim::Int)
    qvs = get_qvs(qcp)
    i = 1
    while index(qvs) != dim
        i += 1
    end
    return quantifier(qcp.qvs[i])
end