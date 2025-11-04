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


n = 1                                                   # 1
p = 5                                                   # Dimension of the domain of f
P_x(t) = sin(t+π/4) + 2*sin(3*t-3*π/4-1) + sin(3.2*t+π/4-0.9)
P_y(t) = cos(t+π/4) + 2*cos(3*t-3*π/4-1) + cos(3.2*t+π/4-0.9)
f_fun = [x -> sqrt((x[1] - P_x(x[3]))^2 + (x[2] - P_y(x[3]))^2) - x[4] - x[5]]               # function f
Df_fun = [x -> [(x[1] - sin(π/4 + x[3]) - 2sin(-3*π/4 -1 + 3x[3]) - sin(π/4 -0.9 + 3.2x[3])) / sqrt((x[1] - sin(π/4 + x[3]) - 2sin(-3*π/4 -1 + 3x[3]) - sin(π/4 -0.9 + 3.2x[3]))^2 + (x[2] + cos(π/4 -0.9 + 3.2x[3]) - cos(π/4 + x[3]) + 2cos(-3*π/4 -1 + 3x[3]))^2) (x[2] + cos(π/4 -0.9 + 3.2x[3]) - cos(π/4 + x[3]) + 2cos(-3*π/4 -1 + 3x[3])) / sqrt((x[1] - sin(π/4 + x[3]) - 2sin(-3*π/4 -1 + 3x[3]) - sin(π/4 -0.9 + 3.2x[3]))^2 + (x[2] + cos(π/4 -0.9 + 3.2x[3]) - cos(π/4 + x[3]) + 2cos(-3*π/4 -1 + 3x[3]))^2) (2(x[1] - sin(π/4 + x[3]) - 2sin(-3*π/4 -1 + 3x[3]) - sin(π/4 -0.9 + 3.2x[3]))*(-3.2cos(π/4 -0.9 + 3.2x[3]) - cos(π/4 + x[3]) - 6cos(-3*π/4 -1 + 3x[3])) + 2(x[2] + cos(π/4 -0.9 + 3.2x[3]) - cos(π/4 + x[3]) + 2cos(-3*π/4 -1 + 3x[3]))*(sin(π/4 + x[3]) - 6sin(-3*π/4 -1 + 3x[3]) - 3.2sin(π/4 -0.9 + 3.2x[3]))) / (2sqrt((x[1] - sin(π/4 + x[3]) - 2sin(-3*π/4 -1 + 3x[3]) - sin(π/4 -0.9 + 3.2x[3]))^2 + (x[2] + cos(π/4 -0.9 + 3.2x[3]) - cos(π/4 + x[3]) + 2cos(-3*π/4 -1 + 3x[3]))^2)) -1 -1]]
problem = Problem(f_fun, Df_fun)                        # problem := f, Df
qe = QuantifierProblem(problem, [(Forall, 3), (Exists, 4), (Exists, 5)], p, n)    # problem, [∃ z], p, n
X_0 = IntervalBox(interval(0, 5), interval(0, 5))     # Domain to be paved X_0
intervals = [interval(0, 2), interval(1,20), interval(0, 0)]                           # z ∈ Z = [0, 16]

using TimerOutputs
const to = TimerOutput()

eps = 0.1                                               # Precision
@timeit to "pave not refined" inn, out, delta = pave_11(qe, X_0, intervals, eps, is_refined=false)
@timeit to "pave refined" inn, out, delta = pave_11(qe, X_0, intervals, eps, is_refined=true)

print_inn_out_delta(inn, out, delta)

# Show the timing results
time_without = TimerOutputs.time(to["pave not refined"])
time_with = TimerOutputs.time(to["pave refined"])
speedup = time_without / time_with
println(round(time_without / 10^9, digits=3), " seconds without refinement")
println(round(time_with / 10^9, digits=3), " seconds with refinement")
println(round(speedup, digits=2), " speedup due to refinement")

# Display the results on a plot if requested
if parsed_args["display"]
    display(X_0, inn, out, delta)
end