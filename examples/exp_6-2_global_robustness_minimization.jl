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
        "--save"
            help = "save the output"
            action = :store_true
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
# ------------------------------------------------------
ϵ_x = [0.01]
ϵ_p = [0.01, 1]
allow_exists_and_forall_bisection = parsed_args["refine"]
allow_exists_or_forall_bisection = parsed_args["subdivide"]


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
    return [x[1] + x[2], x_h, x_r]
end

function gradient_perturbation(x)
    return [1 1;
            0 0;
            0 0]
end

# Set the problem

p = 3
f_fun = Function[(x...) -> (confidence(N([x[2], x_h, x_r])) - x[1]),
            (x...) -> (confidence(N(perturbation(x[2:3]))))]
# Df_fun = RowJacobianMatrix([(x...) -> [-1 confidence(DN(IntervalBox([x[2], x_h, x_r])) * gradient_variables)... 0],
#             (x...) -> [0 confidence(DN(IntervalBox(perturbation(x[2:3]))) * gradient_perturbation(x[2:3]))...]])
Df_fun = ColJacobianMatrix([(x...) -> [-1, 0],
            (x...) -> [confidence(dN(1, [x[2], x_h, x_r])), confidence(dN(1, perturbation(x[2:3])))],
            (x...) -> [confidence(dN(1, [x[2], x_h, x_r])), confidence(dN(1, perturbation(x[2:3])))]])
sizes = [1, 1]
formula = parseformula("(¬ P_1) ∨ P_2")
qvs = [(Forall, 2), (Forall, 3)]
qvs_relaxed =  [qvs, qvs]
parameters = ProblemParameters(formula, f_fun, Df_fun, sizes, qvs, qvs_relaxed, p)
X_0 = IntervalBox(interval(0.95, 1))
P = [interval(25/149, 1), interval(-2/149, 2/149)]
G = [interval(minus_inf, -strict_epsilon), interval(strict_epsilon, plus_inf)]
domains = ProblemDomains(P, G)

configuration = PavingConfiguration(ϵ_x, ϵ_p, allow_exists_or_forall_bisection, allow_exists_and_forall_bisection)

# Pave

verbose = false
benchmark = @benchmarkable (global inn, out, delta = pave_12(X_0, parameters, domains, configuration, verbose)) samples=1 evals=1
results = run(benchmark)
timing = round(median(results).time/1e6, digits=3) # time/1e6 -> timing in milliseconds

# Save the paving in .png file

filename = splitext(PROGRAM_FILE)[1]
if allow_exists_and_forall_bisection
    base = "$(filename)_$(ϵ_x)_$(ϵ_p)_refined"
elseif allow_exists_or_forall_bisection
    base = "$(filename)_$(ϵ_x)_$(ϵ_p)_subdivided"
else
    base = "$(filename)_$(ϵ_x)"
end

outfile = base * ".png"
save_drawing(X_0, inn, out, delta, outfile)

# Save the timing in a log file

logfile = base * ".log"
save_log(timing, logfile)
