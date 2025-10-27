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


n = 1                                                                           # 1
p = 4                                                                           # Dimension of the domain of f
@variables x[1:p]                                                               # x[1] := x, x[2] := y, x[3] := z, x[4] := ζ
f_num = [x[1]^2+x[2]^2+2*x[1]*x[2]-20*x[1]-20x[2]+100-x[3]-x[4]]                # f(x, y, z, zeta) = x^2 + y^2 + 2xy - 20x - 20y + 100 - z - ζ
f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
problem = Problem(f_fun, Df_fun)                                                # problem := f, Df
qe = QuantifierProblem(problem, [(Forall, 3), (Exists, 2), (Exists, 4)], p, n)  # problem, [∀ z, ∃ y, ∃ ζ], p, n
X_0 = IntervalBox(interval(0, 0.1))                                             # Domain to be paved X_0
intervals = [interval(2, 8), interval(6, 8), interval(0, 0)]                    # z ∈ [2, 8], y ∈ [-0.1, 0.1], ζ ∈ [0, 0]

using TimerOutputs
const to = TimerOutput()

eps = 0.05                                                                      # Precision
# @timeit to "pave not refined" inn, out, delta = pave_11(qe, X_0, intervals, eps, is_refined=false)
@timeit to "pave refined" inn, out, delta = pave_11(qe, X_0, intervals, eps, 1.0, 1, is_refined=true)

print_inn_out_delta(inn, out, delta)

# Show the timing results
show(to)
println()

# Display the results on a plot if requested
if parsed_args["display"]
    display(X_0, inn, out, delta)
end