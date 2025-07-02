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

function negation(qv::QuantifiedVariable)
    quantifier, index = qv
    @match quantifier begin
        $Forall => return (Exists, index)
        $Exists => return (Forall, index)
        _ => error("Unknown quantifier: $quantifier")
    end
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
    quantifiers::Vector{QuantifiedVariable}
    p::Int
    n::Int
end

function QuantifierProblem(f, Df, quantifiers::Vector{Any}, p::Int, n::Int)
    quantifier_variables = QuantifiedVariable[]
    for i in 1:2:length(quantifiers)
        if quantifiers[i] == "forall" && quantifiers[i+1] isa Int
            push!(quantifier_variables, (Forall, quantifiers[i+1]))
        elseif quantifiers[i] == "exists" && quantifiers[i+1] isa Int
            push!(quantifier_variables, (Exists, quantifiers[i+1]))
        elseif quantifiers[i] != "forall" && quantifiers[i] != "exists"
            error("""Invalid quantifier: quantifiers[$i], "$(quantifiers[i])", should be "forall" or "exists".""")
        else
            error("""Invalid quantifier: quantifier[$(i+1)], "$(quantifiers[i+1])", should be an integer.""")
        end
    end
    return QuantifierProblem(f, Df, quantifier_variables, p, n)
end

function quantifier(qe::QuantifierProblem, dim::Int)
    i = 1
    while index(qe.quantifiers[i]) != dim
        i += 1
    end
    return quantifier(qe.quantifiers[i])
end