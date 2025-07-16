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

struct QuantifierProblem
    f::Vector{Function}
    Df::Vector{Function}
    qvs::Vector{QuantifiedVariable}
    p::Int
    n::Int
end

function QuantifierProblem(f, Df, qvs::Vector{Any}, p::Int, n::Int)
    @assert length(qvs) % 2 == 0 "Quantifier variables should be in pairs of (quantifier, index)."
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
    return QuantifierProblem(f, Df, quantifier_variables, p, n)
end

function quantifier(qe::QuantifierProblem, dim::Int)
    i = 1
    while index(qe.qvs[i]) != dim
        i += 1
    end
    return quantifier(qe.qvs[i])
end