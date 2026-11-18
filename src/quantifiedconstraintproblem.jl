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

function _build_ranges(sizes::Vector{Int})
    end_indices::Vector{Int} = cumsum(sizes)
    start_indices::Vector{Int} = end_indices .- sizes .+ 1
    return range.(start_indices, end_indices)
end

function build_ranges_dict(M::SyntaxTree, sizes::Vector{Int})
    keys = map(atom -> atom.value, atoms(M))
    ranges = _build_ranges(sizes)
    return Dict(zip(keys, ranges))
end

function build_positions_dict(M::SyntaxTree)
    atom_from_position = Dict(enumerate(atoms(M)))
    position_from_atom = Dict(value => key for (key, value) in atom_from_position)
    return position_from_atom
end

struct ProblemParameters
    M::SyntaxTree
    f::Union{UserDefinedFunctions, UserDefinedSymbolicFunctions}
    sizes::Vector{Int}
    qvs::Vector{QuantifiedVariable}
    p::Int
    function ProblemParameters(M::SyntaxTree, variables::AbstractVector{Num}, f_num::Vector{Num}, sizes::Vector{Int}, qvs::Vector{QuantifiedVariable}, p::Int)
        @assert length(sizes) == length(atoms(M)) "Number of dimensions should be equal to number of predicates."
        userfunctions = UserDefinedSymbolicFunctions(variables, f_num)
        new(M, userfunctions, sizes, qvs, p)
    end
    function ProblemParameters(M::SyntaxTree, variables::AbstractVector{Num}, f_num::Vector{Num}, sizes::Vector{Int}, qvs::Vector{Any}, p::Int)
        ProblemParameters(M, variables, f_num, sizes, convert(Vector{QuantifiedVariable}, qvs), p)
    end
    function ProblemParameters(variables::AbstractVector{Num}, f_num::Vector{Num}, qvs::Vector{QuantifiedVariable}, p::Int)
        n = length(f_num)
        M = parseformula("P")
        sizes = [n]
        ProblemParameters(M, variables, f_num, sizes, qvs, p)
    end
    function ProblemParameters(variables::AbstractVector{Num}, f_num::Vector{Num}, qvs::Vector{Any}, p::Int)
        ProblemParameters(variables, f_num, convert(Vector{QuantifiedVariable}, qvs), p)
    end
    function ProblemParameters(M::SyntaxTree, f::Vector{Function}, Df::JacobianMatrix, sizes::Vector{Int}, qvs::Vector{QuantifiedVariable}, p::Int)
        @assert length(sizes) == length(atoms(M)) "Number of dimensions should be equal to number of predicates."
        userfunctions = UserDefinedFunctions(f, Df)
        new(M, userfunctions, sizes, qvs, p)
    end
    function ProblemParameters(M::SyntaxTree, f::Vector{Function}, Df::JacobianMatrix, sizes::Vector{Int}, qvs::Vector{Any}, p::Int)
        ProblemParameters(M, f, Df, sizes, convert(Vector{QuantifiedVariable}, qvs), p)
    end
end

get_M(parameters::ProblemParameters) = parameters.M
get_f(parameters::ProblemParameters) = parameters.f
get_sizes(parameters::ProblemParameters) = parameters.sizes
get_qvs(parameters::ProblemParameters) = parameters.qvs
get_p(parameters::ProblemParameters) = parameters.p

struct ProblemDomains{T}
    P::Vector{IntervalArithmetic.Interval{T}}
    G::Vector{IntervalArithmetic.Interval{T}}
    function ProblemDomains(P, G)
        new{Float64}(convert(Vector{IntervalArithmetic.Interval{Float64}}, P), convert(Vector{IntervalArithmetic.Interval{Float64}}, G))
    end
end

get_P(domains::ProblemDomains) = domains.P
get_G(domains::ProblemDomains) = domains.G

struct ExpandedProblem{T}
    M::SyntaxTree
    f::Vector{Function}
    Df::JacobianMatrix
    G::Vector{IntervalArithmetic.Interval{T}}
    positions_dict::Dict{Atom{String}, Int}
    function ExpandedProblem(M::SyntaxTree, f::Vector{Function}, Df::JacobianMatrix, G::Vector{IntervalArithmetic.Interval{T}}, positions_dict::Dict{Atom{String}, Int}) where T <: Number
        @assert allunique(atoms(M)) "All predicates in the formula must be unique."
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

struct UnderapproxQuantifiedConstraintProblem
    problem::ExpandedProblem
    qvs::Vector{QuantifiedVariable}
    qvs_relaxed::Vector{Vector{QuantifiedVariable}}
    function UnderapproxQuantifiedConstraintProblem(problem::ExpandedProblem, qvs::Vector{QuantifiedVariable}, qvs_relaxed::Vector{Vector{QuantifiedVariable}})
        @assert length(get_f(problem)) == length(qvs_relaxed) "Number of functions in the problem should be equal to the number of relaxed quantifier variable lists."
        new(problem, qvs, qvs_relaxed)
    end
    function UnderapproxQuantifiedConstraintProblem(problem::ExpandedProblem, qvs::Vector{Any}, qvs_relaxed::Vector{Vector{Any}})
        UnderapproxQuantifiedConstraintProblem(problem, convert(Vector{QuantifiedVariable}, qvs), convert(Vector{Vector{QuantifiedVariable}}, qvs_relaxed))
    end
end

@forward((UnderapproxQuantifiedConstraintProblem, :problem), ExpandedProblem)
get_qvs(qcp::UnderapproxQuantifiedConstraintProblem) = qcp.qvs
get_qvs_relaxed(qcp::UnderapproxQuantifiedConstraintProblem) = qcp.qvs_relaxed

struct OverapproxQuantifiedConstraintProblem
    problem::ExpandedProblem
    qvs::Vector{QuantifiedVariable}
    function OverapproxQuantifiedConstraintProblem(problem::ExpandedProblem, qvs::Vector{QuantifiedVariable})
        new(problem, qvs)
    end
    function OverapproxQuantifiedConstraintProblem(problem::ExpandedProblem, qvs::Vector{Any})
        OverapproxQuantifiedConstraintProblem(problem, convert(Vector{QuantifiedVariable}, qvs))
    end
end

@forward((OverapproxQuantifiedConstraintProblem, :problem), ExpandedProblem)
get_qvs(qcp::OverapproxQuantifiedConstraintProblem) = qcp.qvs
