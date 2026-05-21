<<<<<<<< HEAD:examples/ex_5-5_robot_collision.jl
include("../src/pave.jl")
========
include("../src/utils.jl")
include("../src/pave.jl")
include("../src/read_onnx.jl")
>>>>>>>> cav:examples/exp_6-2_global_robustness_lambda_delta_optimization.jl

using BenchmarkTools

# ------------------------------------------------------
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
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
ϵ_x = [0.1, 0.1]
ϵ_p = [0.01, 1]
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

nnet = read_onnx_mlp("normed_etcs_2x250_nodes.onnx")
N(x) = NeuralVerification.compute_output(nnet, x)
DN(x) = get_gradient(nnet, x)

confidence(vect::AbstractVector) = vect[1] - vect[2]
confidence(mat::Matrix) = mat[1, :] - mat[2, :]

gradient_variables = [1; 0; 0]

x_h = 20/50
x_r = 35/50

function perturbation(x)
    return [x[3] + x[4] * x[1], x_h, x_r]
end

function gradient_perturbation(x)
    return [x[4] 0 1 x[1]; 0 0 0 0; 0 0 0 0]
end

n = 2
p = 6
f_fun = [x -> (confidence(N([x[3], x_h, x_r])) - x[2] - x[5]), x -> (confidence(N(perturbation(x[1:4]))) - x[6])]
Df_fun = [x -> [0 -1 confidence(DN(IntervalBox([x[3], x_h, x_r])) * gradient_variables)... 0 -1 0],
        x -> [confidence(DN(IntervalBox(perturbation(x[1:4]))) * gradient_perturbation(x[1:4]))... 0 -1]]
problem = Problem(f_fun, Df_fun, [[1], [2]])
qvs = [(Forall, 3), (Forall, 4)]
qcp = QuantifiedConstraintProblem(problem, qvs, [qvs, qvs], p, n)
X_0 = IntervalBox(interval(0, 1), interval(0, 1))
p_in = [[interval(100/149, 1)], [interval(-0.013, 0.013)]]
p_out = deepcopy(p_in)
G = [interval(minus_inf, -strict_epsilon), interval(strict_epsilon, plus_inf)]

if allow_exists_and_forall_bisection
    refine_in!(p_in, qcp.qvs, ϵ_p, qcp.p, qcp.n)
    refine_out!(p_out, qcp.qvs, ϵ_p, qcp.p, qcp.n)
end

using TimerOutputs
const to = TimerOutput()

<<<<<<<< HEAD:examples/ex_5-5_robot_collision.jl
@btime (global inn, out, delta = pave_11(X_0, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection)) samples=10
println("Undecided domain: ", round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")
========
@timeit to "pave" inn, out, delta = pave_monotonous_12(X_0, [1, -1], p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection)

show(to)
>>>>>>>> cav:examples/exp_6-2_global_robustness_lambda_delta_optimization.jl

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
<<<<<<<< HEAD:examples/ex_5-5_robot_collision.jl
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
========

    if allow_exists_and_forall_bisection
        outfile = "$(filename)_$(ϵ_x)_$(ϵ_p)_refined.png"
    elseif allow_exists_or_forall_bisection
        outfile = "$(filename)_$(ϵ_x)_$(ϵ_p)_subdivided.png"
    else
        outfile = "$(filename)_$(ϵ_x).png"
    end
    savefig(outfile)
    println("The result was saved in $(outfile).")
end
>>>>>>>> cav:examples/exp_6-2_global_robustness_lambda_delta_optimization.jl
