include("../src/utils.jl")
include("../src/pave.jl")

using NeuralVerification
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
        "--refine", "-r"
            help = "bisect the parameters with ∀ and replace the ones with ∃ by points (or vice versa), requires eps_p, does not work with --subdivide"
            action = :store_true
        "--subdivide", "-s"
            help = "bisect the parameters with either ∀ or ∃, requires eps_p, does not work with --refine"
            action = :store_true
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
# ------------------------------------------------------

filename = splitext(PROGRAM_FILE)[1]

# Load the neural network

nnet = read_nnet("models/running_example.nnet"; last_layer_activation = NeuralVerification.ReLU())
N(x) = NeuralVerification.compute_output(nnet, x)
DN(x) = get_gradient(nnet, x)

# Define the confidence function

confidence(vect::AbstractVector) = vect[1] - vect[2]
confidence(mat::Matrix) = mat[1, :] - mat[2, :]

# Set the problem

p = 2
f_fun = Function[x -> confidence(N(x[1:2]))]
Df_fun = Function[x -> [confidence(DN(x[1:2]))...]]
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
ϵ_p = nothing
allow_exists_and_forall_bisection = parsed_args["refine"]
allow_exists_or_forall_bisection = parsed_args["subdivide"]
configuration = PavingConfiguration(ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection)

# Pave

inn, out, delta = pave_12(X_0, parameters, domains, configuration)

# Save the paving in .png file

if get_allow_exists_and_forall_bisection(configuration)
    outfile = "$(filename)_$(δ)_$(ϵ_x)_$(ϵ_p)_refined.png"
elseif get_allow_exists_or_forall_bisection(configuration)
    outfile = "$(filename)_$(δ)_$(ϵ_x)_$(ϵ_p)_subdivided.png"
else
    outfile = "$(filename)_$(δ)_$(ϵ_x).png"
end

save_drawing(X_0, inn, out, delta, outfile)