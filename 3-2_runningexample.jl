#!/usr/local/bin/julia

# using IntervalArithmetic

include("genreach2.jl")
include("pave.jl")
include("quantifierproblem.jl")

# ------------------------------------------------------
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--display", "-d"
            help = "display the results in a plot (to be used in interactive mode with 'julia -i ex_runningexample.jl -d')"
            action = :store_true
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
# ------------------------------------------------------


n = 1                                                               # 1
p = 3                                                               # Dimension of the domain of f
@variables x[1:p]                                                   # x[1] := x, x[2] := p_1, x[3] := z
f_num = [x[2]^2-(x[1]-1)*(x[1]-2)*(x[1]-3)-x[3]]                    # f(x, p_1, z) = p_1^2 - (x - 1)(x - 2)(x - 3) - z
f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
problem = Problem(f_fun, Df_fun)                                    # problem := f, Df
qe = QuantifierProblem(problem, [(Forall, 2), (Exists, 3)], p, n)   # problem, [∀ p1, ∃ z], p, n
X_0 = IntervalBox(interval(0, 5))                                   # Domain to be paved X_0
intervals = [interval(0, 1/4), interval(-1/4, 1/4)]                 # p_1 ∈ [0, 1/4], z ∈ [-1/4, 1/4]

using TimerOutputs
const to = TimerOutput()

eps = 0.1                                                                # Precision
@timeit to "pave not refined" inn, out, delta = pave_11(qe, X_0, intervals, eps, is_refined=false)
@timeit to "pave refined" inn, out, delta = pave_11(qe, X_0, intervals, eps, is_refined=true)

print_inn_out_delta(inn, out, delta)

# Show the timing results
show(to)
println()

# Display the results on a plot if requested
if parsed_args["display"]
    display(X_0, inn, out, delta)
end