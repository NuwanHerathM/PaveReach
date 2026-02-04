include("pave.jl")

# ------------------------------------------------------
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
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
f_num = [x[1]^2 + x[2]^2 - x[3]]
f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
problem = Problem(f_fun, Df_fun)
qvs = []
qcp = QuantifiedConstraintProblem(problem, qvs, [qvs], p, n)
X_0 = IntervalBox(interval(-5, 5), interval(-5, 5))
p_in = []
p_out = deepcopy(p_in)
G = [interval(0, 16)]

ϵ_x = 1
ϵ_p = 0.5
allow_exists_or_forall_bisection = false
allow_exists_and_forall_bisection = false

@assert nand(allow_exists_and_forall_bisection, allow_exists_or_forall_bisection)

if allow_exists_and_forall_bisection
    refine_in!(p_in, qcp.qvs, ϵ_p, qcp.p, qcp.n)
    refine_out!(p_out, qcp.qvs, ϵ_p, qcp.p, qcp.n)
end

println("ϵ_x  = ", ϵ_x)
println("ϵ_p  = ", ϵ_p)
println(if allow_exists_and_forall_bisection "Refined" else "Not refined" end)
println(if allow_exists_or_forall_bisection "Normal bisection on P" else "No standard bisection on P" end)

inn, out, delta = pave_11(X_0, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection)
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
    width = 1000
    height = 1000
    buffer = 50
    Drawing(width + 2*buffer, height + 2*buffer, outfile)
    luxor_draw(X_0, inn, out, delta, width, height, buffer)
    finish()
    println("The result was saved in $(outfile).")
end