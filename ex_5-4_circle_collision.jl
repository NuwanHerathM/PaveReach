include("pave.jl")

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
        "eps_p"
            help = "epsilon for the parameters"
            arg_type = Float64
            required = false
        "--refine", "-r"
            help = "bisect the parameters with ∀ and replace the ones with ∃ by points (or vice versa), requires eps_p, does not work with --subdivide"
            action = :store_true
        "--subdivide", "-s"
            help = "bisect the parameters with either ∀ or ∃, requires eps_p, does not work with --refine"
            action = :store_true
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
ϵ_x = parsed_args["eps_x"]
ϵ_p = parsed_args["eps_p"]
allow_exists_and_forall_bisection = parsed_args["refine"]
allow_exists_or_forall_bisection = parsed_args["subdivide"]


if isnothing(ϵ_p)
    @assert !allow_exists_and_forall_bisection "eps_p was not provided. Refinement requires eps_p. Use --help for more information."
    @assert !allow_exists_or_forall_bisection "eps_p was not provided. Subdivision requires eps_p. Use --help for more information."
else
    if !allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection
        @warn "eps_p was provided, but will not be used. Provide options --refine or --subdivide in order to use eps_p."
    end
end
@assert nand(allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) "Refinement and subdivision are mutually exclusive. Use --help for more information."

use_plots = parsed_args["with_plots"]
use_luxor = parsed_args["with_luxor"]

@assert nand(use_plots, use_luxor) "Plots and Luxor are mutually exclusive. Use --help for more information."
# ------------------------------------------------------

filename = splitext(PROGRAM_FILE)[1]

n = 1
p = 4
f_fun = [x -> (2.5*sin(x[3]) - x[1])^2 + (2.5*cos(x[3]) - x[2])^2 - x[4]]
Df_fun = [x -> [(-2*(2.5*sin(x[3]) - x[1])) (-2*(2.5*cos(x[3]) - x[2])) (5*cos(x[3])*(2.5*sin(x[3]) - x[1]) - 5*sin(x[3])*(2.5*cos(x[3]) - x[2])) -1]]
problem = Problem(f_fun, Df_fun)
qvs = [(Forall, 3)]
qcp = QuantifiedConstraintProblem(problem, qvs, [qvs], p, n)
X_0 = IntervalBox(interval(-5, 5), interval(-5, 5))
p_in = [[interval(-π, π)]]
p_out = deepcopy(p_in)
G = [interval(0.25, 100)]

if allow_exists_and_forall_bisection
    refine_in!(p_in, qcp.qvs, ϵ_p, qcp.p, qcp.n)
    refine_out!(p_out, qcp.qvs, ϵ_p, qcp.p, qcp.n)
end

println("ϵ_x = ", ϵ_x)
println("ϵ_p = ", ϵ_p)

@btime (global inn, out, delta = pave_11(X_0, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection))
println("Undecided domain: ", round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")

if allow_exists_and_forall_bisection
    outfile = "$(filename)_11_$(ϵ_x)_$(ϵ_p)_refined.png"
elseif allow_exists_or_forall_bisection
    outfile = "$(filename)_11_$(ϵ_x)_$(ϵ_p)_subdivided.png"
else
    outfile = "$(filename)_11_$(ϵ_x).png"
end

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