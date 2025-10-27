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


n = 1                                                        # 1
p = 3                                                        # Dimension of the domain of f
@variables x[1:p]                                            # x[1] := x_1, x[2] := x_2, x[3] := z
f_num_1 = [x[1]^2 + x[2]^2 - x[3]]                           # f_1(x, z) = x_1^2 + x_2^2 - z
f_fun_1, Df_fun_1 = build_function_f_Df(f_num_1, x, n, p)
problem_1 = Problem(f_fun_1, Df_fun_1)                       # problem_1 := f_1, Df_1
f_num_2 = [(2 - x[1])^2 + x[2]^2 - x[3]]                     # f_2(x, z) = x_1^2 + x_2^2 - z
f_fun_2, Df_fun_2 = build_function_f_Df(f_num_2, x, n, p)
problem_2 = Problem(f_fun_2, Df_fun_2)                       # problem_2 := f_2, Df_2
problem = AndProblem([problem_1, problem_2])                 # problem = problem_1 ∧ problem_2
qe = QuantifierProblem(problem, [(Exists, 3)], p, n)         # problem, [∃ z], p, n
X_0 = IntervalBox(interval(-5, 5), interval(-5, 5))          # Domain to be paved X_0
intervals = [interval(0, 16)]                                # z ∈ Z = [0, 16]

using TimerOutputs
const to = TimerOutput()

eps = 0.1                                                    # Precision
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