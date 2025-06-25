#!/usr/local/bin/julia

using Match

@enum Quantifier begin
    Forall
    Exists
end

QuantifiedVariable = Tuple{Quantifier, Int}

function quantifier(qv::QuantifiedVariable)::Quantifier
    return qv[1]
end

function index(qv::QuantifiedVariable)::Int
    return qv[2]
end

function print(qv::QuantifiedVariable)
    quantifier, index = qv
    @match quantifier begin
        $Forall => Base.print("∀", index)
        $Exists => Base.print("∃", index)
        _ => error("Unknown quantifier: $quantifier")
    end
end

function println(qv::QuantifiedVariable)
    print(qv)
    Base.println()
end

function print(qvs::Vector{QuantifiedVariable})
    Base.print("[")
    for (i, qv) in enumerate(qvs)
        print(qv)
        if i < length(qvs)
            Base.print(", ")
        end
    end
    Base.print("]")
end

function println(qvs::Vector{QuantifiedVariable})
    print(qvs)
    Base.println()
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
    fun::Vector{Num}
    quantifiers::Vector{QuantifiedVariable}
    p::Int
    n::Int
end

function QuantifierProblem(fun::Vector{Num}, quantifiers::Vector{Any}, p::Int, n::Int)
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
    return QuantifierProblem(fun, quantifier_variables, p, n)
end

function quantifier(qe::QuantifierProblem, dim::Int)::Quantifier
    i = 1
    while index(qe.quantifiers[i]) != dim
        i += 1
    end
    return quantifier(qe.quantifiers[i])
end