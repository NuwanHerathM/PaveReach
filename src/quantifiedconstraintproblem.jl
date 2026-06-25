#!/usr/local/bin/julia

using Match
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

function get_negated_predicates(M::SyntaxTree)
    l = SyntaxTree[M]
    negated = String[]
    while !isempty(l)
        tree = pop!(l)
        if token(tree) == ¬
            if length(tree.children) == 1 && first(tree.children) isa Atom
                push!(negated, first(tree.children).value)
            else
                error("Negation of non-atomic formula is not supported.")
            end
        elseif (token(tree) == ∧) || (token(tree) == ∨)
            append!(l, tree.children)
        elseif tree isa Atom
            continue
        else
            error("Unknown syntax tree tree type: $(typeof(tree))")
        end        
    end
    return negated
end

function get_indices_range(sizes)
    end_indices::Vector{Int} = cumsum(sizes)
    start_indices::Vector{Int} = end_indices .- sizes .+ 1
    return range.(start_indices, end_indices)
end

function get_indices_dict(M, sizes)
    keys = map(atom -> atom.value, atoms(M))
    ranges = get_indices_range(sizes)
    return Dict(zip(keys, ranges))
end

function get_positions_dict(M)
    atom_from_position = Dict(enumerate(atoms(M)))
    position_from_atom = Dict(value => key for (key, value) in atom_from_position)
    return position_from_atom
end

struct Problem
    M::SyntaxTree
    f::Vector{Function}
    Df::Vector{Function}
    positions_dict::Dict{Atom{String}, Int}
    function Problem(M::SyntaxTree, f, Df, positions_dict::Dict{Atom{String}, Int})
        @assert length(f) == length(Df) "Number of functions and gradients must be equal."
        @assert length(f) == maximum(values(positions_dict)) "Number of functions and largest position must be equal."
        @assert allunique(atoms(M)) "All predicates in the formula must be unique."
        @assert length(positions_dict) == natoms(M) "Number of dictionary entries must match the number of predicates in the formula."
        new(M, f, Df, positions_dict)
    end
    # function Problem(formula_string::String, f::Vector{Function}, Df::Vector{Function}, sizes::Vector{Int})
    #     @assert length(f) == length(Df) "Number of functions and gradients must be equal."
    #     @assert length(f) == sum(sizes) "Number of functions and dimensions must be equal."
    #     M = parseformula(formula_string)
    #     @assert allunique(atoms(M)) "All predicates in the formula must be unique."
    #     @assert length(sizes) == natoms(M) "Number of dimensions must match the number of predicates in the formula."
    #     indices_dict = get_indices_dict(M, sizes)
    #     negated_predicates = get_negated_predicates(M)
    #     new(M, f, Df, indices_dict, negated_predicates)
    # end
end

# abstract type ConnectedProblem end

# struct AndProblem <: ConnectedProblem
#     problems::Vector{Union{Problem, ConnectedProblem}}
# end

# struct OrProblem <: ConnectedProblem
#     problems::Vector{Union{Problem, ConnectedProblem}}
# end

# problems(problem::ConnectedProblem) = problem.problems
# problems(problem::Problem) = [problem]

struct QuantifiedConstraintProblem
    problem::Problem
    qvs::Vector{QuantifiedVariable}
    qvs_relaxed::Vector{Vector{QuantifiedVariable}}
    p::Int
    n::Int
    function QuantifiedConstraintProblem(problem::Problem, qvs, qvs_relaxed, p::Int, n::Int)
        @assert length(problem.f) == n "Number of functions in the problem should be equal to n = $n."
        @assert length(qvs_relaxed) == n "Number of relaxed quantifier variable lists should be equal to n = $n."
        new(problem, qvs, qvs_relaxed, p, n)
    end
end

function QuantifiedConstraintProblem(f, Df, dnf_indices, qvs::Vector{Any}, p::Int, n::Int)
    @assert length(qvs) % 2 == 0 "Quantifier variables should be in pairs of (quantifier, index)."
    problem = Problem(f, Df, dnf_indices)
    for i in 1:2:length(qvs)
        quantifier = qvs[i]
        idx = qvs[i+1]
        if quantifier == "forall" && idx isa Int
            push!(quantifier_variables, (Forall, idx))
        elseif quantifier == "exists" && idx isa Int
            push!(quantifier_variables, (Exists, idx))
        elseif quantifier != "forall" && quantifier != "exists"
            error("""Invalid quantifier: qvs[$i], "$(quantifier)", should be "forall" or "exists".""")
        else
            error("""Invalid quantifier: qvs[$(i+1)], "$(idx)", should be an integer.""")
        end
    end
    return QuantifiedConstraintProblem(problem, quantifier_variables, [quantifier_variables], p, n)
end

function quantifier(qcp::QuantifiedConstraintProblem, dim::Int)
    i = 1
    while index(qcp.qvs[i]) != dim
        i += 1
    end
    return quantifier(qcp.qvs[i])
end