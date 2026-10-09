include("../src/PaveReach.jl")

using BenchmarkTools

# Load the neural network

nnet = read_nnet("models/running_example.nnet"; last_layer_activation = NeuralVerification.ReLU())
N(x) = NeuralVerification.compute_output(nnet, x)
DN(x) = get_gradient(nnet, x)
dN(i, x) = get_partial_derivative(nnet, i, x)
N(x::Tuple) = N(collect(x))
DN(x::Tuple) = DN(collect(x))
dN(i, x::Tuple) = dN(i, collect(x))


# Define the confidence function

confidence(vect::AbstractVector) = vect[1] - vect[2]
confidence(mat::Matrix) = mat[1, :] - mat[2, :]

function perturbation(x)
    return [x[1] + x[3], x[2] + x[4]]
end

function gradient_perturbation(x)
    return [1 0 1 0;
            0 1 0 1]
end

# Set the problem

p = 4
f_fun = Function[(x...) -> confidence(N(x[1:2])),
            (x...) -> confidence(N(perturbation(x[1:4])))]
df_1(x...) = [confidence(dN(1, x[1:2])), confidence(dN(1, perturbation(x[1:4])))]
df_2(x...) = [confidence(dN(2, x[1:2])), confidence(dN(2, perturbation(x[1:4])))]
df_3(x...) = [0, confidence(dN(1, perturbation(x[1:4])))]
df_4(x...) = [0, confidence(dN(2, perturbation(x[1:4])))]
Df_fun = ColJacobianMatrix([df_1, df_2, df_3, df_4])
sizes = [1, 1]
formula = parseformula("(¬ P_1) ∨ P_2")
qvs = [(Forall, 3), (Forall, 4)]
qvs_relaxed = [qvs, qvs]
parameters = ProblemParameters(formula, f_fun, Df_fun, sizes, qvs, qvs_relaxed, p)
X_0 = IntervalBox(interval(-1, 1), interval(-1, 1))
ϵ_max = 1/8
P = [interval(-ϵ_max, ϵ_max), interval(-ϵ_max, ϵ_max)]
δ = 1/8
G = [interval(δ + strict_epsilon, plus_inf), interval(strict_epsilon, plus_inf)]
domains = ProblemDomains(P, G)

ϵ_x = [0.05, 0.05]
configuration = PavingConfiguration(ϵ_x)

# Pave

verbose = false
benchmark = @benchmarkable (global inn, out, delta = pave_12(X_0, parameters, domains, configuration, verbose)) samples=1 evals=1
results = run(benchmark)
timing = round(median(results).time/1e6, digits=3) # time/1e6 -> timing in milliseconds

# Save the paving in .png file

filename = splitext(PROGRAM_FILE)[1]
base = "$(filename)_$(δ)_$(ϵ_x)"
    
outfile = base * ".png"
save_drawing(X_0, inn, out, delta, outfile)

# Save the timing in a log file

logfile = base * ".log"
save_log(timing, logfile)