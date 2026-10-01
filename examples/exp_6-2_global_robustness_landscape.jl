include("../src/utils.jl")
include("../src/pave.jl")
include("../src/read_onnx.jl")

using BenchmarkTools

# ------------------------------------------------------
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "delta"
            help = "confidence delta"
            arg_type = Float64
            required = true
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
# ------------------------------------------------------

# Load the neural network

nnet = read_onnx_mlp("models/normed_etcs_2x250_nodes.onnx")
N(x) = NeuralVerification.compute_output(nnet, x)
DN(x) = get_gradient(nnet, x)

x_h = 20/50
x_r = 35/50

# Define the confidence function

confidence(vect::AbstractVector) = vect[1] - vect[2]
confidence(mat::Matrix) = mat[1, :] - mat[2, :]

gradient_variables = [1; 0; 0]

function perturbation(x)
    return [x[1] + x[2], x_h, x_r]
end

function gradient_perturbation(x)
    return [1 1; 0 0; 0 0]
end

# Set the problem

p = 2
f_fun = [x -> confidence(N([x[1], x_h, x_r])), x -> confidence(N(perturbation(x[1:2])))]
Df_fun = [x -> [confidence(DN(IntervalBox([x[1], x_h, x_r])) * gradient_variables)... 0],
        x -> [confidence(DN(IntervalBox(perturbation(x[1:2]))) * gradient_perturbation(x[1:2]))...]]
sizes = [1, 1]
formula = parseformula("(¬ P_1) ∨ P_2")
qvs = [(Forall, 2)]
qvs_relaxed = [qvs, qvs]
parameters = ProblemParameters(formula, f_fun, Df_fun, sizes, qvs, qvs_relaxed, p)
X_0 = IntervalBox(interval(25/149, 1))
P = [interval(-2/149, 2/149)]
δ = parsed_args["delta"]
G = [interval(δ, plus_inf), interval(strict_epsilon, plus_inf)]
domains = ProblemDomains(P, G)

ϵ_x = [0.01]
configuration = PavingConfiguration(ϵ_x)

# Pave

# @btime(global inn, out, delta = pave_12(X_0, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection))
benchmark = @benchmarkable (global inn, out, delta = pave_12(X_0, parameters, domains, configuration)) samples=1 evals=1
results = run(benchmark)
timing = round(median(results).time/1e6, digits=3)
@info "Time (ms)" t=timing

# Save the paving in .png file

filename = splitext(PROGRAM_FILE)[1]
outfile = "$(filename)_$(δ)_$(ϵ_x).png"

xticks = (25, 149)
save_drawing(X_0, inn, out, delta, outfile, xticks)