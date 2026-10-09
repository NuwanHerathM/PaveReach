include("../src/PaveReach.jl")

using BenchmarkTools

# Load the neural network

nnet = read_onnx_mlp("models/normed_etcs_2x250_nodes.onnx")
N(x) = NeuralVerification.compute_output(nnet, x)
DN(x) = get_jacobian(nnet, x)
dN(i, x) = get_partial_derivative(nnet, i, x)
N(x::Tuple) = N(collect(x))
DN(x::Tuple) = DN(collect(x))
dN(i, x::Tuple) = dN(i, collect(x))

# Define the confidence function

confidence(vect::AbstractVector) = vect[2] - vect[1]
confidence(mat::Matrix) = mat[2, :] - mat[1, :]

# Set the problem

p = 3
f_fun = Function[(x...) -> confidence(N([x[1], x[2], x[3]]))]
{{jacobian}}
sizes = [1]
formula = parseformula("P_1")
qvs = [(Forall, 2), (Forall, 3)]
qvs_relaxed = [qvs]
parameters = ProblemParameters(formula, f_fun, Df_fun, sizes, qvs, qvs_relaxed, p)
X_0 = IntervalBox(interval(25/149, 1))
P = [interval(0, 800/50000), interval(0, 800/50000)]
δ = 0.0001
G = [interval(δ, plus_inf)]
domains = ProblemDomains(P, G)

ϵ_x = [0.01]
configuration = PavingConfiguration(ϵ_x)

# Pave

verbose = false
benchmark = @benchmarkable (global inn, out, delta = pave_12(X_0, parameters, domains, configuration, verbose)) samples=1 evals=1
results = run(benchmark)
timing = round(median(results).time/1e6, digits=3) # time/1e6 -> timing in milliseconds

# Save the paving in .png file

filename = splitext(PROGRAM_FILE)[1]
base = "$(filename)_$(δ)_$(ϵ_x)"

xticks = (25, 149)
outfile = base * ".png"
save_drawing(X_0, inn, out, delta, outfile, xticks)

# Save the timing in a log file

logfile = base * ".log"
save_log(timing, logfile)