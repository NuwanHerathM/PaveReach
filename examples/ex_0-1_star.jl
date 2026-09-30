include("../src/pave.jl")

filename = splitext(PROGRAM_FILE)[1]

# Set the problem

p = 2
# x[1] := x, x[2] := y
@variables x[1:p]
# f_num := f
f_num = [(x[1] + 1)^2 + x[2]^2,                                                                 # f_1(x, y) = (x + 1)^2 + y^2
    (x[1] - 1)^2 + x[2]^2,                                                                      # f_2(x, y) = (x - 1)^2 + y^2
    x[1]^2 + (x[2] + 1)^2,                                                                      # f_3(x, y) = x^2 + (y + 1)^2
    x[1]^2 + (x[2] - 1)^2]                                                                      # f_4(x, y) = x^2 + (y - 1)^2
sizes = [2, 2]                                                                                  # P_1 <- (f_1, f_2), P_2 <- (f_3, f_4)
formula = parseformula("P_1 ∨ P_2")
qvs = QuantifiedVariable[]                                                                      # No quantified variables in this example
qvs_relaxed = [qvs, qvs, qvs, qvs]
# 
parameters = ProblemParameters(formula, x, f_num, sizes, qvs, qvs_relaxed, p)
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

inn, out, delta = pave(X_0, parameters, domains, configuration)
println("Undecided domain: ", round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")

# Save the paving in .png file

outfile = "$(filename)_$(ϵ_x).png"

save_drawing(X_0, inn, out, delta, outfile)