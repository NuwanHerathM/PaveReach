include("../src/PaveReach.jl")

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

# Load the neural network

nnet = read_onnx_mlp("models/normed_etcs_2x250_nodes.onnx")
N(x) = NeuralVerification.compute_output(nnet, x)
DN(x) = get_jacobian(nnet, x)
dN(i, x) = get_partial_derivative(nnet, i, x)
N(x::Tuple) = N(collect(x))
DN(x::Tuple) = DN(collect(x))
dN(i, x::Tuple) = dN(i, collect(x))

# Define the confidence function

confidence(vect::AbstractVector) = vect[1] - vect[2]
confidence(mat::Matrix) = mat[1, :] - mat[2, :]

gradient_variables = [1; 0; 0]

x_h = 20/50
x_r = 35/50

function perturbation(x)
    return [x[3] + x[4] * x[1], x_h, x_r]
end

function gradient_perturbation(x)
    return [x[4] 0 1 x[1];
            0 0 0 0;
            0 0 0 0]
end

# Set the problem
# { λ, δ ∣ ∀ v, ∀ ϵ, confidence(N(v, x_h, x_r)) > δ ⟹ confidence(N(v + λϵ, x_h, x_r)) > 0 }

p = 4
# x[1] = λ, x[2] = δ, x[3] = v, x[4] = ϵ
f_fun = Function[(x...) -> confidence(N([x[3], x_h, x_r])) - x[2],
            (x...) -> confidence(N(perturbation(x[1:4])))]
Df_fun = RowJacobianMatrix([(x...) -> [0, -1, confidence(DN(IntervalBox([x[3], x_h, x_r])) * gradient_variables)..., 0],
            (x...) -> [confidence(DN(IntervalBox(perturbation(x[1:4]))) * gradient_perturbation(x[1:4]))...]])
Df_fun = ColJacobianMatrix([(x...) -> [0, confidence(dN())]])
sizes = [1, 1]
formula = parseformula("(¬ P_1) ∨ P_2")
qvs = [(Forall, 3), (Forall, 4)]
qvs_relaxed = [qvs, qvs]
parameters = ProblemParameters(formula, f_fun, Df_fun, sizes, qvs, qvs_relaxed, p)
X_0 = IntervalBox(interval(0, 1), interval(0, 1))
P = [interval(100/149, 1), interval(-0.013, 0.013)]
G = [interval(strict_epsilon, plus_inf), interval(strict_epsilon, plus_inf)]
domains = ProblemDomains(P, G)

ϵ_x = [0.1, 0.1]
ϵ_p = [0.01, 1]
allow_exists_and_forall_bisection = parsed_args["refine"]
allow_exists_or_forall_bisection = parsed_args["subdivide"]
configuration  = PavingConfiguration(ϵ_x, ϵ_p, allow_exists_or_forall_bisection, allow_exists_and_forall_bisection)

# Pave

inn, out, delta = pave_11(X_0, parameters, domains, configuration)
println("Undecided domain: ", round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")
# results = run(benchmark)
# timing = round(median(results).time/1e6, digits=3) # time/1e6 -> timing in milliseconds

# Save the paving in .png file

filename = splitext(PROGRAM_FILE)[1]
if allow_exists_and_forall_bisection
    base = "$(filename)_11_$(ϵ_x)_$(ϵ_p)_refined"
elseif allow_exists_or_forall_bisection
    base = "$(filename)_11_$(ϵ_x)_$(ϵ_p)_subdivided"
else
    base = "$(filename)_11_$(ϵ_x)"
end

outfile = base * ".png"
save_drawing(X_0, inn, out, delta, outfile)

# Save the timing in a log file

logfile = base * ".log"
save_log(timing, logfile)