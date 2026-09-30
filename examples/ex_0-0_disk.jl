include("../src/pave.jl")

# Set the problem

p = 2
@variables x[1:p]
f_num = [x[1]^2 + x[2]^2]
qvs = []
qvs_relaxed = [qvs]
parameters = ProblemParameters(x, f_num, qvs, qvs_relaxed, p)
X_0 = IntervalBox(interval(-5, 5), interval(-5, 5))
P = []
G = [interval(0, 16)]
domains = ProblemDomains(P, G)

ϵ_x = [0.1, 0.1]
configuration = PavingConfiguration(ϵ_x)

println(configuration)

# Pave

inn, out, delta = pave_11(X_0, parameters, domains, configuration)
println("Undecided domain: ", round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")

# Save the paving in .png file

filename = splitext(PROGRAM_FILE)[1]
outfile = "$(filename)_11_$(ϵ_x).png"

save_drawing(X_0, inn, out, delta, outfile)