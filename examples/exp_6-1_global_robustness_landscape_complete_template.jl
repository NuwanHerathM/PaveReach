include("../src/PaveReach.jl")

using BenchmarkTools

# Load the neural network

nnet = read_nnet("models/running_example.nnet"; last_layer_activation = NeuralVerification.ReLU())
N(x) = NeuralVerification.compute_output(nnet, x)
DN(x) = get_jacobian(nnet, x)
dN(x, i) = get_partial_derivative(nnet, i, x)

# Define the confidence function

confidence_1(vect::AbstractVector) = vect[1] - vect[2]
confidence_1(mat::Matrix) = mat[1, :] - mat[2, :]
confidence_2(vect::AbstractVector) = vect[2] - vect[1]
confidence_2(mat::Matrix) = mat[2, :] - mat[1, :]

function perturbation(x)
    return [x[1] + x[3], x[2] + x[4]]
end

function gradient_perturbation(x)
    return [1 0 1 0; 0 1 0 1]
end

# Set the problem

p = 4
f_fun = Function[x -> confidence_1(N(x[1:2])),
        x -> confidence_1(N(perturbation(x[1:4]))),
        x -> confidence_2(N(x[1:2])),
        x -> confidence_2(N(perturbation(x[1:4])))]
Df_fun = RowJacobianMatrix([x -> [confidence_1(DN(x[1:2]))..., 0, 0],
        x -> [confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4]))...],
        x -> [confidence_2(DN(x[1:2]))..., 0, 0],
        x -> [confidence_2(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4]))...]])
df_1(x) = [confidence_1(dN(x[1:2], 1)), confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 1]), confidence_1(dN(x[1:2], 1)), confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 1])]
df_2(x) = [confidence_1(dN(x[1:2], 2)), confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 2]), confidence_1(dN(x[1:2], 2)), confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 2])]
df_3(x) = [0, confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 3]), 0, confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 3])]
df_4(x) = [0, confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 4]), 0, confidence_1(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4])[:, 4])]
# df_1(x) = [confidence_1(dN(x[1:2], 1)), confidence_1(dN(perturbation(x[1:4]), 1)), confidence_1(dN(x[1:2], 1)), confidence_1(dN(perturbation(x[1:4]), 1))]
# df_2(x) = [confidence_1(dN(x[1:2], 2)), confidence_1(dN(perturbation(x[1:4]), 2)), confidence_1(dN(x[1:2], 2)), confidence_1(dN(perturbation(x[1:4]), 2))]
# df_3(x) = [0, confidence_1(dN(perturbation(x[1:4]), 1)), 0, confidence_1(dN(perturbation(x[1:4]), 1))]
# df_4(x) = [0, confidence_1(dN(perturbation(x[1:4]), 2)), 0, confidence_1(dN(perturbation(x[1:4]), 2))]
Df_fun = ColJacobianMatrix([df_1, df_2, df_3, df_4])
sizes = [1, 1, 1, 1]
formula = parseformula("((¬ P_1) ∨ P_2) ∧ ((¬ P_3) ∨ P_4)")
qvs = [(Forall, 3), (Forall, 4)]
qvs_relaxed = [qvs, qvs, qvs, qvs]
parameters = ProblemParameters(formula, f_fun, Df_fun, sizes, qvs, [qvs, qvs, qvs, qvs], p)
X_0 = IntervalBox(interval(-1, 1), interval(-1, 1))
ϵ_max = 1/8
P = [interval(-ϵ_max, ϵ_max), interval(-ϵ_max, ϵ_max)]
δ = 1/8
G = [interval(δ + strict_epsilon, plus_inf), interval(strict_epsilon, plus_inf), interval(δ + strict_epsilon, plus_inf), interval(strict_epsilon, plus_inf)]
domains = ProblemDomains(P, G)

ϵ_x = [0.05, 0.05]
configuration = PavingConfiguration(ϵ_x)

println(configuration)

# Pave

using TimerOutputs
const to = TimerOutput()

benchmark = @benchmarkable (global inn, out, delta = pave_12(X_0, parameters, domains, configuration)) samples=1 evals=1
results = run(benchmark)
timing = round(median(results).time/1e6, digits=3)
@info "Time (ms)" t=timing

# Save the paving in .png file

filename = splitext(PROGRAM_FILE)[1]
outfile = "$(filename)_$(δ)_$(ϵ_x).png"

save_drawing(X_0, inn, out, delta, outfile)
