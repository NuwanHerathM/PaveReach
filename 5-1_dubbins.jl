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


n = 1                                                                                    # 1
p = 14                                                                                   # Dimension of the domain of f
@variables x[1:p]                                                                        # x[1] := x, x[2] := y, x[3] := θ, x[4:10] := ϵ, x[11] := t, x[12:15] := z
f_num_1 = [0.1*x[4]+(1+0.01*x[5])*x[11]+1.31*10^(-7)*x[6]*x[11]^2-x[1]-x[12]]            # f(x, p_1, z) = x_1^2 + b*x_2^2 - z
f_fun_1, Df_fun_1 = build_function_f_Df(f_num_1, x, n, p)
problem_1 = Problem(f_fun_1, Df_fun_1)                                                   # problem := f, Df
f_num_2 = [0.1*x[7]+(0.01*x[9]+0.01*x[10]*x[11])*x[11]+(0.005*x[8])*x[11]^2-x[2]-x[13]]  # f(x, p_1, z) = x_1^2 + b*x_2^2 - z
f_fun_2, Df_fun_2 = build_function_f_Df(f_num_2, x, n, p)
problem_2 = Problem(f_fun_2, Df_fun_2)                                                   # problem := f, Df
f_num_3 = [0.01*x[9]+0.01*x[10]*x[11]-x[3]-x[14]]                                        # f(x, p_1, z) = x_1^2 + b*x_2^2 - z
f_fun_3, Df_fun_3 = build_function_f_Df(f_num_3, x, n, p)
problem_3 = Problem(f_fun_3, Df_fun_3)                                                   # problem := f, Df
problem = AndProblem([problem_1, problem_2, problem_3])                                  # problem := f, Df
qvs = [(Exists, 10), (Exists, 4), (Exists, 7), (Exists, 9), (Forall, 5), (Forall, 6), (Forall, 8), (Exists, 11), (Exists, 12), (Exists, 13), (Exists, 14)]  # [∃ ϵ_7, ∃ ϵ_1, ∃ ϵ_4, ∃ ϵ_6, ∀ ϵ_2, ∀ ϵ_3, ∀ ϵ_5, ∃ t, ∃ z_1, ∃ z_2, ∃ z_3]
qe = QuantifierProblem(problem, qvs, p, n)                                               # problem, qvs, p, n
X_0 = IntervalBox(interval(-1, 1), interval(-0.5, 0.5), interval(-0.05,0.05))                                      # Domain to be paved X_0
intervals = [repeat([interval(-1, 1)], 7)..., interval(0, 0.5),repeat([interval(-2, 2)], 3)...]                                        # z ∈ Z = [0, 16], b ∈ [-0.1, 0.1]

using TimerOutputs
const to = TimerOutput()

eps = 0.1                                                                # Precision
# @timeit to "pave not refined" inn, out, delta = pave_11(qe, X_0, intervals, eps, is_refined=false)
@timeit to "pave refined" inn, out, delta = pave_11(qe, X_0, intervals, eps, is_refined=true)

print_inn_out_delta(inn, out, delta)

# Show the timing results
show(to)
println()

# Display the results on a plot if requested
if parsed_args["display"]
    display(X_0, inn, out, delta)
end