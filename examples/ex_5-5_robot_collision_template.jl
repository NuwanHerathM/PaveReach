include("../src/PaveReach.jl")

using Plots.PlotMeasures
using LaTeXStrings
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
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
# ------------------------------------------------------

# Set the problem

p = 3
P_x(t) = sin(t+π/4) + 2*sin(3*t-3*π/4-1) + sin(3.2*t+π/4-0.9)
P_y(t) = cos(t+π/4) + 2*cos(3*t-3*π/4-1) + cos(3.2*t+π/4-0.9)
dP_x(t) = cos(t+π/4) + 6*cos(3*t-3*π/4-1) + 3.2*cos(3.2*t+π/4-0.9)
dP_y(t) = -sin(t+π/4) - 6*sin(3*t-3*π/4-1) - 3.2*sin(3.2*t+π/4-0.9)
f_fun = Function[(x...) -> (x[1] - P_x(x[3]))^2 + (x[2] - P_y(x[3]))^2]
{{jacobian}}
sizes = [1]
formula = parseformula("P_1")
qvs = [(Forall, 3)]
qvs_relaxed = [qvs]
parameters = ProblemParameters(formula, f_fun, Df_fun, sizes, qvs, qvs_relaxed, p)
X_0 = IntervalBox(interval(0, 5), interval(0, 5))
P = [interval(0, 2)]
G = [interval(0.25, 25)]
domains = ProblemDomains(P, G)

ϵ_x = [0.1, 0.1]
ϵ_p = [0.1]
allow_exists_and_forall_bisection = true
allow_exists_or_forall_bisection = false
configuration = PavingConfiguration(ϵ_x, ϵ_p, allow_exists_or_forall_bisection, allow_exists_and_forall_bisection)

# Pave

@btime (global inn, out, delta = pave_11(X_0, parameters, domains, configuration)) samples=1 evals=1
println("Undecided domain: ", round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")

# Save the paving in .png file

filename = splitext(PROGRAM_FILE)[1]
if allow_exists_and_forall_bisection
    outfile = "$(filename)_11_$(ϵ_x)_$(ϵ_p)_refined.png"
elseif allow_exists_or_forall_bisection
    outfile = "$(filename)_11_$(ϵ_x)_$(ϵ_p)_subdivided.png"
else
    outfile = "$(filename)_11_$(ϵ_x).png"
end

save_drawing(X_0, out, inn, delta, outfile)