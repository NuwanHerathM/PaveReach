include("../src/pave.jl")

using Plots.PlotMeasures
using LaTeXStrings
using Match
using BenchmarkTools

# ------------------------------------------------------
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "o_in"
            help = "oracle for IN"
            arg_type = Int
            required = true
        "o_out"
            help = "oracle for OUT"
            arg_type = Int
            required = true
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
o_in = parsed_args["o_in"]
o_out = parsed_args["o_out"]

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
p = 3
@variables x[1:p]
f_num = [x[2]^2-(x[1]-1)*(x[1]-2)*(x[1]-3)-x[3]]
f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
problem = Problem(f_fun, Df_fun) 
qcp = QuantifiedConstraintProblem(problem, [(Forall, 2)], [[(Forall, 2)]], p, n)
X_0 = IntervalBox(interval(0, 5))
p_in = [[interval(0, 1/4)]]
p_out = deepcopy(p_in)
G = [interval(-1/4, 1/4)]


if allow_exists_and_forall_bisection
    refine_in!(p_in, qcp.qvs, ϵ_p, qcp.p, qcp.n)
    refine_out!(p_out, qcp.qvs, ϵ_p, qcp.p, qcp.n)
end

println("O^IN_$(o_in), O^OUT_$(o_out)")
println("ϵ_x  = ", ϵ_x)
println("ϵ_p = ", ϵ_p)
println(if allow_exists_and_forall_bisection "Refined" else "Not refined" end)
println(if allow_exists_or_forall_bisection "Standard bisection on P" else "No standard bisection on P" end)


paving_function = @match (o_in, o_out) begin
    (1, 1) => pave_11
    (1, 2) => pave_12
    (2, 1) => pave_21
    (2, 2) => pave_22
end

if o_in == 1
    lemma_in = @match (allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) begin
        (false, false) => L"Lemma~3.1"
        (false, true)  => L"Lemma~4.2"
        (true, false)  => L"Lemma~4.4"
    end
else
    lemma_in = @match (allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) begin
        (false, false) => L"Lemma~3.3"
        (false, true)  => L"a la Lemma~4.2"
        (true, false)  => L"a la Lemma~4.4"
    end
end
if o_out == 1
    lemma_out = @match (allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) begin
        (false, false) => L"Lemma~3.2"
        (false, true)  => L"Lemma~4.3"
        (true, false)  => L"Lemma~4.5"
    end
else
    lemma_out = @match (allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) begin
        (false, false) => L"Lemma~3.4"
        (false, true)  => L"a la Lemma~4.3"
        (true, false)  => L"a la Lemma~4.5"
    end
end

@btime (global inn, out, delta = paving_function(X_0, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection)) samples=10
println("Undecided domain: ", round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")

if allow_exists_and_forall_bisection
    outfile = "$(filename)_$(o_in)$(o_out)_$(ϵ_x)_$(ϵ_p)_refined.png"
elseif allow_exists_or_forall_bisection
    outfile = "$(filename)_$(o_in)$(o_out)_$(ϵ_x)_$(ϵ_p)_subdivided.png"
else
    outfile = "$(filename)_$(o_in)$(o_out)_$(ϵ_x).png"
end

if use_plots
    p = plot()
    plot_title = L"$\mathcal{O}^{IN}$ " * lemma_in * L", $\mathcal{O}^{OUT}$ " * lemma_out
    title!(p, plot_title)
    draw(p, X_0, inn, out, delta)
    plot(p, size=(600, 80), grid=false, titlelocation = :left, topmargin=5Plots.PlotMeasures.mm, titlefontsize=12)
    savefig(outfile)
    println("The result was saved in $(outfile).")
end

if use_luxor
    width = 600
    height = 40
    buffer = 50
    Drawing(width + 2*buffer, height + 2*buffer, outfile)
    luxor_draw(X_0, inn, out, delta, width, height, buffer)
    fontsize(16)
    plot_title = L"\mathcal{O}^{IN}~%$(lemma_in), \mathcal{O}^{OUT}~%$(lemma_out)"
    Luxor.text(plot_title, Point(buffer, buffer/2), valign = :middle)
    finish()
    println("The result was saved in $(outfile).")
end