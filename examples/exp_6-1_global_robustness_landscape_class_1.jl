include("../src/utils.jl")
include("../src/pave.jl")

using BenchmarkTools

# Load the neural network

nnet = read_nnet("models/running_example.nnet"; last_layer_activation = NeuralVerification.ReLU())
N(x) = NeuralVerification.compute_output(nnet, x)
DN(x) = get_gradient(nnet, x)

# Define the confidence function

confidence(vect::AbstractVector) = vect[1] - vect[2]
confidence(mat::Matrix) = mat[1, :] - mat[2, :]

function perturbation(x)
    return [x[1] + x[3], x[2] + x[4]]
end

function gradient_perturbation(x)
    return [1 0 1 0; 0 1 0 1]
end

# Set the problem

p = 4
f_fun = Function[x -> confidence(N(x[1:2])), x -> confidence(N(perturbation(x[1:4])))]
Df_fun = Function[x -> [confidence(DN(x[1:2]))... 0 0],
        x -> [confidence(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4]))...]]
sizes = [1, 1]
formula = parseformula("(¬ P_1) ∨ P_2")
qvs = [(Forall, 3), (Forall, 4)]
qvs_relaxed = [qvs, qvs]
parameters = ProblemParameters(formula, f_fun, Df_fun, sizes, qvs, qvs_relaxed, p)
X_0 = IntervalBox(interval(-1, 1), interval(-1, 1))
ϵ_max = 1/8
P = [interval(-ϵ_max, ϵ_max), interval(-ϵ_max, ϵ_max)]
δ = parsed_args["delta"]
G = [interval(δ + strict_epsilon, plus_inf), interval(strict_epsilon, plus_inf)]
domains = ProblemDomains(P, G)

ϵ_x = [0.05, 0.05]
configuration = PavingConfiguration(ϵ_x)

# Pave

using TimerOutputs
const to = TimerOutput()

@timeit to "pave" inn, out, delta = pave_12(X_0, parameters, domains, configuration)

show(to)
println()

# Save the paving in .png file

filename = splitext(PROGRAM_FILE)[1]
outfile = "$(filename)_$(δ)_$(ϵ_x).png"
    
save_drawing(X_0, inn, out, delta, outfile)