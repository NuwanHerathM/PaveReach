include("../src/PaveReach.jl")

using BenchmarkTools

input_x = 0.5
input_y = -0.35

filename = splitext(PROGRAM_FILE)[1]

nnet = read_nnet("models/running_example.nnet"; last_layer_activation = NeuralVerification.ReLU())
N(x) = NeuralVerification.compute_output(nnet, x)
DN(x) = get_jacobian(nnet, x)
dN(i, x) = get_partial_derivative(nnet, i, x)
N(x::Tuple) = N(collect(x))
DN(x::Tuple) = DN(collect(x))
dN(i, x::Tuple) = dN(i, collect(x))

confidence(vect::AbstractVector) = vect[1] - vect[2]
confidence(mat::Matrix) = mat[1, :] - mat[2, :]

input = [input_x, input_y]

function perturbation(x)
    return [input[1] + x[1] * x[2], input[2] + x[1] * x[3]]
end

function gradient_perturbation(x)
    return [x[2] x[1] 0;
            x[3] 0 x[1]]
end

p = 3
f_fun = Function[(x...) -> confidence(N(perturbation(x[1:3])))]
Df_fun = ColJacobianMatrix([(x...) -> [confidence(DN(perturbation(x[1:3])) * [x[2]; x[3]])],
                                (x...) -> [confidence(dN(1, perturbation(x[1:3])) * x[1])],
                                (x...) -> [confidence(dN(2, perturbation(x[1:3])) * x[1])]])
sizes = [1]
formula = parseformula("P_1")
qvs = [(Forall, 2), (Forall, 3)]
qvs_relaxed = [qvs]
parameters = ProblemParameters(formula, f_fun, Df_fun, sizes, qvs, qvs_relaxed, p)
X_0 = IntervalBox(interval(0, 1))
ϵ_max = 0.5
P = [interval(-ϵ_max, ϵ_max), interval(-ϵ_max, ϵ_max)]
G = [interval(strict_epsilon, 100000)]
domains = ProblemDomains(P, G)

ϵ_x = [0.05]
ϵ_p = [0.1, 0.1]
allow_exists_or_forall_bisection = false
allow_exists_and_forall_bisection = true
configuration = PavingConfiguration(ϵ_x, ϵ_p, allow_exists_or_forall_bisection, allow_exists_and_forall_bisection)

# Pave

verbose = false
benchmark = @benchmarkable (global inn, out, delta = pave_12(X_0, parameters, domains, configuration, verbose)) samples=1 evals=1
results = run(benchmark)
timing = round(median(results).time/1e6, digits=3) # time/1e6 -> timing in milliseconds

# Save the paving in a .png file

base = "$(filename)_$(input_x)_$(input_y)_$(ϵ_x)_$(ϵ_p)_refined"

outfile = base * ".png"
save_drawing(X_0, inn, out, delta, outfile)

# Save the timing in a log file

logfile = base * ".log"
save_log(timing, logfile)