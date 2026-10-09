include("../src/PaveReach.jl")

using BenchmarkTools

filename = splitext(PROGRAM_FILE)[1]

nnet = read_nnet("models/running_example.nnet"; last_layer_activation = NeuralVerification.ReLU())
N(x) = NeuralVerification.compute_output(nnet, x)
DN(x) = get_jacobian(nnet, x)
dN(i, x) = get_partial_derivative(nnet, i, x)
N(x::Tuple) = N(collect(x))
DN(x::Tuple) = DN(collect(x))
dN(i, x::Tuple) = dN(i, collect(x))


confidence_1(vect::AbstractVector) = vect[1] - vect[2]
confidence_1(mat::Matrix) = mat[1, :] - mat[2, :]
confidence_2(vect::AbstractVector) = vect[2] - vect[1]
confidence_2(mat::Matrix) = mat[2, :] - mat[1, :]

function perturbation(x)
    return [x[1] + x[3], x[2] + x[4]]
end

function gradient_perturbation(x)
    return [1 0 1 0;
            0 1 0 1]
end

p = 5
f_fun = [(x...) -> (confidence_1(N(x[2:3])) - x[1]),
            (x...) -> (confidence_1(N(perturbation(x[2:5])))),
            (x...) -> (confidence_2(N(x[2:3])) - x[1] ),
            (x...) -> (confidence_2(N(perturbation(x[2:5]))))]
Df_fun = RowJacobianMatrix([(x...) -> [-1, confidence_1(DN(x[2:3]))..., 0, 0],
            (x...) -> [0, confidence_1(DN(perturbation(x[2:5])) * gradient_perturbation(x[2:5]))...],
            (x...) -> [-1, confidence_2(DN(x[2:3]))..., 0, 0],
            (x...) -> [0, confidence_2(DN(perturbation(x[2:5])) * gradient_perturbation(x[2:5]))...]])
# Df_fun = ColJacobianMatrix([(x...) -> [-1, 0, -1, 0],
# (x...) -> [confidence_1(dN(1, x[2:3])), confidence_1(dN(1, perturbation(x[2:5]))), confidence_2(dN(1, x[2:3])), confidence_2(dN(1, perturbation(x[2:5])))],
# (x...) -> [confidence_1(dN(2, x[2:3])), confidence_1(dN(2, perturbation(x[2:5]))), confidence_2(dN(2, x[2:3])), confidence_2(dN(2, perturbation(x[2:5])))],
# (x...) -> [0, confidence_1(dN(1, perturbation(x[2:5]))), 0, confidence_2(dN(1, perturbation(x[2:5])))],
# (x...) -> [0, confidence_1(dN(2, perturbation(x[2:5]))), 0, confidence_2(dN(2, perturbation(x[2:5])))]])
sizes = [1, 1, 1, 1]
formula = parseformula("((¬ P_1) ∨ P_2) ∧ ((¬ P_3) ∨ P_4)")
qvs = [(Forall, 2), (Forall, 3), (Forall, 4), (Forall, 5)]
qvs_relaxed = [qvs, qvs, qvs, qvs]
parameters = ProblemParameters(formula, f_fun, Df_fun, sizes, qvs, qvs_relaxed, p)
X_0 = IntervalBox(interval(0, 2))
ϵ_max = 1/8
P = [interval(-1, 1), interval(-1, 1), interval(-ϵ_max, ϵ_max), interval(-ϵ_max, ϵ_max)]
G = [interval(strict_epsilon, plus_inf), interval(strict_epsilon, plus_inf), interval(strict_epsilon, plus_inf), interval(strict_epsilon, plus_inf)]
domains = ProblemDomains(P, G)

ϵ_x = [0.1]
ϵ_p = [0.05, 0.05, 0.1, 0.1]
allow_exists_or_forall_bisection = false
allow_exists_and_forall_bisection = true
configuration = PavingConfiguration(ϵ_x, ϵ_p, allow_exists_or_forall_bisection, allow_exists_and_forall_bisection)

# Pave

# inn, out, delta = pave_monotonous_mid_12(X_0, [-1], parameters, domains, configuration)
inn, out, delta = pave_monotonous_sides_optimized_12(X_0, [-1], parameters, domains, configuration)
# inn, out, delta = pave_12(X_0, parameters, domains, configuration)
# results = run(benchmark)
# timing = round(median(results).time/1e6, digits=3) # time/1e6 -> timing in milliseconds

display(out)

# Save the paving in a .png file

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