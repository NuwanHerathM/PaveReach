include("../src/PaveReach.jl")

using BenchmarkTools

filename = splitext(PROGRAM_FILE)[1]

# Set the problem

p = 2
f_fun = Function[(x...) -> (x[1] + 1)^2 + x[2]^2,
    (x...) -> (x[1] - 1)^2 + x[2]^2,
    (x...) -> x[1]^2 + (x[2] + 1)^2,
    (x...) -> x[1]^2 + (x[2] - 1)^2]
Df_fun = MatJacobianMatrix([(x...) -> 2*(x[1] + 1) (x...) -> 2*x[2];
        (x...) -> 2*(x[1] - 1) (x...) -> 2*x[2];
        (x...) -> 2*x[1] (x...) -> 2*(x[2] + 1);
        (x...) -> 2*x[1] (x...) -> 2*(x[2] - 1)])
sizes = [2, 2]
formula = parseformula("P_1 ∨ P_2")
qvs = QuantifiedVariable[]
qvs_relaxed = [qvs, qvs, qvs, qvs]
# 
parameters = ProblemParameters(formula, f_fun, Df_fun, sizes, qvs, qvs_relaxed, p)
X_0 = IntervalBox(interval(-5, 5), interval(-5, 5))                                             # X_0 = [-5, 5] x [-5, 5]
P = IntervalArithmetic.Interval{Float64}[]                                                      # No parameters in this example
G = [interval(minus_inf, 4),                                                                    # f_1(x, y) <= 4 that is f_1(x, y) ∈ (-∞, 4]
    interval(minus_inf, 4),                                                                     # f_2(x, y) <= 4 that is f_2(x, y) ∈ (-∞, 4]
    interval(minus_inf, 4),                                                                     # f_3(x, y) <= 4 that is f_3(x, y) ∈ (-∞, 4]
    interval(minus_inf, 4)]                                                                     # f_4(x, y) <= 4 that is f_4(x, y) ∈ (-∞, 4]
domains = ProblemDomains(P, G)

ϵ_x = [0.1, 0.1]
configuration = PavingConfiguration(ϵ_x)

println(configuration)

# Pave

@btime (global inn, out, delta = pave(X_0, parameters, domains, configuration)) samples=1 evals=1
println("Undecided domain: ", round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")

# Save the paving in .png file

outfile = "$(filename)_$(ϵ_x).png"

save_drawing(X_0, inn, out, delta, outfile)