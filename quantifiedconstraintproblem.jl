#!/usr/local/bin/julia

using Match

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

struct Problem
    f::Vector{Function}
    Df::Vector{Function}
    ϕ::Vector{Function}
    dnf_1_indices::Vector{Vector{Int}}
    dnf_2_indices::Vector{Vector{Int}}
end

abstract type ConnectedProblem end

struct AndProblem <: ConnectedProblem
    problems::Vector{Union{Problem, ConnectedProblem}}
end

struct OrProblem <: ConnectedProblem
    problems::Vector{Union{Problem, ConnectedProblem}}
end

problems(problem::ConnectedProblem) = problem.problems
problems(problem::Problem) = [problem]

struct QuantifiedConstraintProblem
    problem::Union{Problem, ConnectedProblem}
    qvs::Vector{QuantifiedVariable}
    qvs_relaxed::Vector{Vector{QuantifiedVariable}}
    p::Int
    n::Int
end

function QuantifiedConstraintProblem(f, Df, qvs::Vector{Any}, p::Int, n::Int)
    @assert length(qvs) % 2 == 0 "Quantifier variables should be in pairs of (quantifier, index)."
    problem = Problem(f, Df)
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