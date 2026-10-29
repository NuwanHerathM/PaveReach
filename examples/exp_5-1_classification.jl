include("../src/PaveReach.jl")

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

modelpath = joinpath(@__DIR__, "models/running_example.nnet")
nnet = read_nnet(modelpath; last_layer_activation = NeuralVerification.ReLU())
N(x) = NeuralVerification.compute_output(nnet, x)
DN(x) = get_gradient(nnet, x)
dN(i, x) = get_partial_derivative(nnet, i, x)
N(x::Tuple) = N(collect(x))
DN(x::Tuple) = DN(collect(x))
dN(i, x::Tuple) = dN(i, collect(x))

# Define the confidence function

confidence(vect::AbstractVector) = vect[1] - vect[2]
confidence(mat::Matrix) = mat[1, :] - mat[2, :]

# Set the problem
# x_1, x_2 ∣ confidence(N(x_1, x_2)) > δ

p = 2
f_fun = Function[(x...) -> confidence(N(x[1:2]))]
Df_fun = ColJacobianMatrix([(x...) -> [confidence(dN(i, x[1:2]))] for i in 1:p])
sizes = [1]
formula = parseformula("P_1")
qvs = QuantifiedVariable[]
qvs_relaxed = [qvs]
parameters = ProblemParameters(formula, f_fun, Df_fun, sizes, qvs, qvs_relaxed, p)
X_0 = IntervalBox(interval(-1, 1), interval(-1, 1))
P = []
δ = parsed_args["delta"]
G = [interval(δ + strict_epsilon, plus_inf)]
domains = ProblemDomains(P, G)

ϵ_x = [0.05, 0.05]
configuration = PavingConfiguration(ϵ_x)

# Pave

verbose = false
benchmark = @benchmarkable (global inn, out, delta = pave_12(X_0, parameters, domains, configuration, verbose)) samples=1 evals=1
results = run(benchmark)
timing = round(median(results).time/1e6, digits=3) # time/1e6 -> timing in milliseconds

# Save the paving in .png file

filename = splitext(basename(PROGRAM_FILE))[1]
base = "$(filename)_$(δ)_$(ϵ_x)"

outfile = base * ".png"
save_drawing(X_0, inn, out, delta, @__DIR__, outfile)

# Save the timing in a log file

logfile = base * ".log"
save_log(timing, @__DIR__, logfile)