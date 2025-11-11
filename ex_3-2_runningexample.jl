include("pave2.jl")

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
        "eps"
            help = "epsilon for the paving"
            arg_type = Float64
            required = true
        "factor"
            help = "precision factor on the parameter"
            arg_type = Float64
            required = true
        "--refine", "-r"
            help = "refine the parameters"
            action = :store_true
        "--subdivide", "-s"
            help = "subdivide the parameters"
            action = :store_true
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
# ------------------------------------------------------

filename = splitext(PROGRAM_FILE)[1]

n = 1
p = 3
@variables x[1:p]
f_num = [x[2]^2-(x[1]-1)*(x[1]-2)*(x[1]-3)-x[3]]
f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
problem = Problem(f_fun, Df_fun) 
qe = QuantifierProblem(problem, [(Forall, 2)], [[(Forall, 2)]], p, n)
X_0 = IntervalBox(interval(0, 5))
p_in = [[interval(0, 1/4)]]
p_out = deepcopy(p_in)
G = [interval(-1/4, 1/4)]

o_in = parsed_args["o_in"]
o_out = parsed_args["o_out"]

ϵ = parsed_args["eps"]
factor = parsed_args["factor"]
is_refined = parsed_args["refine"]
allow_normal_p_bisect = parsed_args["subdivide"]
if is_refined
    refine_in!(p_in, qe.qvs, ϵ, qe.p, qe.n)
    refine_out!(p_out, qe.qvs, ϵ, qe.p, qe.n)
end

println("O^IN_$(o_in), O^OUT_$(o_out)")
println("ϵ  = ", ϵ)
println("factor: ", factor)
println(if is_refined "Refined" else "Not refined" end)
println(if allow_normal_p_bisect "Normal bisection on P" else "No standard bisection on P" end)

@assert nand(is_refined, allow_normal_p_bisect)

paving_function = @match (o_in, o_out) begin
    (1, 1) => pave_11
    (1, 2) => pave_12
    (2, 1) => pave_21
    (2, 2) => pave_22
end
plot_title = @match (o_in, o_out) begin
    (1, 1) => L"$\mathcal{O}^{IN}(\mathbb{X}, \mathbb{P}, \mathbb{G}), \mathcal{O}^{OUT}(\mathbb{X}, \mathbb{P}, \mathbb{G})$"
    (1, 2) => L"$\mathcal{O}^{IN}(\mathbb{X}, \mathbb{P}, \mathbb{G}), \mathcal{O}^{OUT}(\mathbb{X}, \neg\mathbb{P}, \mathbb{G}^\complement)$"
    (2, 1) => L"$\mathcal{O}^{IN}(\mathbb{X}, \neg\mathbb{P}, \mathbb{G}^\complement), \mathcal{O}^{OUT}(\mathbb{X}, \mathbb{P}, \mathbb{G})$"
    (2, 2) => L"$\mathcal{O}^{IN}(\mathbb{X}, \neg\mathbb{P}, \mathbb{G}^\complement), \mathcal{O}^{OUT}(\mathbb{X}, \neg\mathbb{P}, \mathbb{G}^\complement)$"
end

@btime (global inn, out, delta = paving_function(X_0, p_in, p_out, G, qe, ϵ, factor, is_refined, allow_normal_p_bisect))
p = plot()
title!(p, plot_title)
draw(p, X_0, inn, out, delta)
println(inn)
# print_inn_out_delta(inn, out, delta)
# print_delta_width(delta)
println(round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")

plot(p, size=(600, 80), grid=false, titlelocation = :left, topmargin=5mm, titlefontsize=12)

str_refined = is_refined ? "_refined" : ""
str_subdivided = allow_normal_p_bisect ? "_subdivided" : ""

savefig("$(filename)_$(o_in)$(o_out)_$(ϵ)$(str_refined)$(str_subdivided).png")