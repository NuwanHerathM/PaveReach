include("../src/pave.jl")

using Plots.PlotMeasures
using LaTeXStrings
using BenchmarkTools

# ------------------------------------------------------
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "eps_x"
            help = "epsilon for the paving"
            arg_type = Float64
            required = true
        "--with_plots"
            help = "generate the output with Plots.jl"
            action = :store_true
        "--with_luxor"
            help = "generate the output with Luxor.jl"
            action = :store_true
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
# ------------------------------------------------------
use_plots = parsed_args["with_plots"]
use_luxor = parsed_args["with_luxor"]

@assert nand(use_plots, use_luxor) "Plots and Luxor are mutually exclusive. Use --help for more information."
# ------------------------------------------------------

filename = splitext(PROGRAM_FILE)[1]

n = 1
p = 3
@variables x[1:p]
# x[1] := v, x[2] := w, x[3] := z
f_num = [-5*x[1]^2-13*x[1]+x[1]*x[2]-x[2]-x[3]]
f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
problem = Problem(f_fun, Df_fun)
qvs = []
qcp = QuantifiedConstraintProblem(problem, qvs, [qvs], p, n)
X_0 = IntervalBox(interval(2, 10), interval(40, 50))
p_in = []
p_out = deepcopy(p_in)
G = [interval(10^(-7), 100)]

ϵ_x = parsed_args["eps_x"]
println("ϵ_x = ", ϵ_x)
ϵ_p = 0.0
is_refined = false
allow_normal_p_bisect = false

@btime (global inn, out, delta = pave_11(X_0, p_in, p_out, G, qcp, ϵ_x, ϵ_p, is_refined, allow_normal_p_bisect)) samples=10
println("Undecided domain: ", round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")

outfile = "$(filename)_$(ϵ_x).png"

if use_plots
    p = plot()
    draw(p, X_0, inn, out, delta)
    plot(p)
    savefig(outfile)
    println("The result was saved in $(outfile).")
end

if use_luxor
    width = 600
    height = 600
    buffer = 50
    Drawing(width + 2*buffer, height + 2*buffer, outfile)
    luxor_draw(X_0, inn, out, delta, width, height, buffer)
    finish()
    println("The result was saved in $(outfile).")
end