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
        "--save", "-s"
            help = "save the plot to a file"
            action = :store_true
        "--time", "-t"
            help = "display the timings"
            action = :store_true
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
# ------------------------------------------------------


n = 1                                                    # 1
p = 3                                                    # Dimension of the domain of f
@variables x[1:p]                                        # x[1] := x_1, x[2] := x_2, x[3] := z
f_num = [x[1]^2 + x[2]^2 - x[3]]                         # f(x, p_1, z) = x_1^2 + x_2^2 - z
f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
problem = Problem(f_fun, Df_fun)                         # problem := f, Df
qvs = [(Exists, 3)]                                     # [∃ z]
qe = QuantifierProblem(problem, qvs, [qvs], p, n)     # problem, [∃ z], p, n
X_0 = IntervalBox(interval(-1, 1), interval(-1, 1))      # Domain to be paved X_0
intervals = [interval(0, 16)]                            # z ∈ Z = [0, 16]

using TimerOutputs
const to = TimerOutput()

filename = splitext(PROGRAM_FILE)[1]

eps = 0.1                                                # Precision
println("ϵ = ", eps)

# @timeit to "pave not refined 1" inn, out, delta = pave_11(qe, X_0, intervals, eps, is_refined=false)
# println(round(volume(delta)/volume(X_0)*100, digits=1), " %")
# draw(X_0, inn, out, delta)
# if parsed_args["save"]
#     savefig("$(filename)_1.png")
# end

@timeit to "pave refined 1" inn, out, delta = pave_11(qe, X_0, intervals, eps, is_refined=true)
println(round(volume(delta)/volume(X_0)*100, digits=1), " %")
draw(X_0, inn, out, delta)
if parsed_args["save"]
    savefig("$(filename)_1_refined.png")
end

show(to)
println()

# if parsed_args["time"]
#     time_without_1 = TimerOutputs.time(to["pave not refined 1"])
#     time_with_1 = TimerOutputs.time(to["pave refined 1"])
#     speedup_1 = time_without_1 / time_with_1
#     println(round(time_without_1 / 10^9, digits=3), " seconds without refinement")
#     println(round(time_with_1 / 10^9, digits=3), " seconds with refinement")
#     println(round(speedup_1, digits=2), " speedup due to refinement")
# end


# @timeit to "pave not refined 2" inn, out, delta = pave_22(qe, X_0, intervals, eps, is_refined=false)
# println(round(volume(delta)/volume(X_0)*100, digits=1), " %")
# draw(X_0, inn, out, delta)
# if parsed_args["save"]
#     savefig("$(filename)_2.png")
# end

# @timeit to "pave refined 2" inn, out, delta = pave_22(qe, X_0, intervals, eps, is_refined=true)
# println(round(volume(delta)/volume(X_0)*100, digits=1), " %")
# draw(X_0, inn, out, delta)
# if parsed_args["save"]
#     savefig("$(filename)_2_refined.png")
# end

# if parsed_args["time"]
#     time_without_2 = TimerOutputs.time(to["pave not refined 2"])
#     time_with_2 = TimerOutputs.time(to["pave refined 2"])
#     speedup_2 = time_without_2 / time_with_2
#     println(round(time_without_2 / 10^9, digits=3), " seconds without refinement")
#     println(round(time_with_2 / 10^9, digits=3), " seconds with refinement")
#     println(round(speedup_2, digits=2), " speedup due to refinement")
# end

# print_inn_out_delta(inn, out, delta)

# Display the results on a plot if requested
if parsed_args["display"]
    gui()
end