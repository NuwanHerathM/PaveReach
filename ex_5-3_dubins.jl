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
        "--x", "-x"
            help = "show problem for x"
            action = :store_true
        "--y", "-y"
            help = "show problem for y"
            action = :store_true
        "--theta", "-t"
            help = "show problem for θ"
            action = :store_true
        "--refine", "-r"
            help = "refine the parameters"
            action = :store_true
        "--subdivide", "-s"
            help = "subdivide the parameters"
            action = :store_true
        "--save"
            help = "save the output"
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
# ------------------------------------------------------

filename = splitext(PROGRAM_FILE)[1]

if parsed_args["x"]
    # x
    n = 1
    p = 6
    # x[1] := x,  x[2:4] := ϵ_1, ϵ_2, ϵ_3, x[5] := t, x[6] := z
    @variables x[1:p]
    f_num_1 = [x[1]-0.1*x[2]-(1+0.01*x[3])*x[5]-1.31*10^(-7)*x[4]*x[5]^2-x[6]]
    f_fun_1, Df_fun_1 = build_function_f_Df(f_num_1, x, n, p)
    problem_1 = Problem(f_fun_1, Df_fun_1)
    qvs_1 = [(Exists, 2), (Forall, 3), (Exists, 4), (Exists, 5)]
    qcp = QuantifiedConstraintProblem(problem_1, qvs_1, [qvs_1], p, n)
    X_0 = IntervalBox(interval(-2, 4))
    p_in = [[interval(-10, 10)], [interval(-1, 1)], [interval(-1, 1)], [interval(0, 2)]]
    p_out = deepcopy(p_in)
    G = [interval(0, 0)]
end

if parsed_args["y"]
    # y
    n = 1 
    p = 7
    # x[1] := y, x[2:5] := ϵ_4, ϵ_5, ϵ_6, ϵ_7, x[6] := t, x[7] := z
    @variables x[1:p]
    f_num_2 = [x[1]-0.1*x[2]-(0.01*x[4]+0.01*x[5]*x[6])*x[6]-(0.005*x[3])*x[6]^2-x[7]]
    f_fun_2, Df_fun_2 = build_function_f_Df(f_num_2, x, n, p)
    problem_2 = Problem(f_fun_2, Df_fun_2)
    qvs_2 = [(Exists, 6), (Exists, 2), (Exists, 4), (Exists, 3), (Exists, 6)]
    qcp = QuantifiedConstraintProblem(problem_2, qvs_2, [qvs_2], p, n)
    X_0 = IntervalBox(interval(-4, 4))
    p_in = [[interval(-10, 10)], [interval(-1, 1)], [interval(-20, 20)], [interval(-1, 1)], [interval(0, 2)]]
    p_out = deepcopy(p_in)
    G = [interval(0, 0)]
end

if parsed_args["theta"]
    # θ
    n = 1
    p = 5
    # x[1] := θ, x[2:3] := ϵ_6, ϵ_7, x[4] := t, x[5] := z
    @variables x[1:p]
    f_num_3 = [x[1]-0.01*x[2]-0.01*x[3]*x[4]-x[5]]
    f_fun_3, Df_fun_3 = build_function_f_Df(f_num_3, x, n, p)
    problem_3 = Problem(f_fun_3, Df_fun_3)
    qvs_3 = [(Exists, 3), (Exists, 2), (Exists, 4)]
    qcp = QuantifiedConstraintProblem(problem_3, qvs_3, [qvs_3], p, n)
    X_0 = IntervalBox(interval(-0.5, 0.5))
    p_in = [[interval(-20,20)], [interval(-1, 1)], [interval(0, 2)]] 
    p_out = deepcopy(p_in)
    G = [interval(0, 0)]
end

if allow_exists_and_forall_bisection
    refine_in!(p_in, qcp.qvs, ϵ_p, qcp.p, qcp.n)
    refine_out!(p_out, qcp.qvs, ϵ_p, qcp.p, qcp.n)
end

println("ϵ_x = ", ϵ_x)
println("ϵ_p = ", ϵ_p)

@btime (global inn, out, delta = pave_11(X_0, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection     ))
print_inn_out_delta(inn, out, delta)
println("Undecided domain: ", round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")

if parsed_args["save"]
    p = plot()
    draw(p, X_0, inn, out, delta)

    if parsed_args["x"]
        title!(p, L"$x$")
        plot(p, size=(600, 80), grid=false, titlelocation = :left, topmargin=5mm, titlefontsize=14)
        outfile = "$(filename)_$(ϵ_x)_$(ϵ_p)_x.png"
    end
    if parsed_args["y"]
        title!(p, L"$y$")
        plot(p, size=(600, 80), grid=false, titlelocation = :left, topmargin=5mm, titlefontsize=14)
        outfile = "$(filename)_$(ϵ_x)_$(ϵ_p)_y.png"
    end
    if parsed_args["theta"]
        title!(p, L"$\theta$")
        plot(p, size=(600, 80), grid=false, titlelocation = :left, topmargin=5mm, titlefontsize=14)
        outfile = "$(filename)_$(ϵ_x)_$(ϵ_p)_theta.png"
    end
    savefig(outfile)
end