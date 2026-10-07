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
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
# ------------------------------------------------------

# Set the problem

p = 3
f_fun = Function[(x...) -> (2.5*sin(x[3]) - x[1])^2 + (2.5*cos(x[3]) - x[2])^2]
{{jacobian}}
sizes = [1]
formula = parseformula("P_1")
qvs = [(Forall, 3)]
qvs_relaxed = [qvs]
parameters = ProblemParameters(formula, f_fun, Df_fun, sizes, qvs, qvs_relaxed, p)
X_0 = IntervalBox(interval(-5, 5), interval(-5, 5))
P = [interval(-π, π)]
G = [interval(0.25, 100)]
domains = ProblemDomains(P, G)

ϵ_x = [0.1, 0.1]
ϵ_p = [0.1]
allow_exists_and_forall_bisection = true
allow_exists_or_forall_bisection = false
configuration = PavingConfiguration(ϵ_x, ϵ_p, allow_exists_or_forall_bisection, allow_exists_and_forall_bisection)

println(configuration)

# Pave

@btime (global inn, out, delta = pave_11(X_0, parameters, domains, configuration)) samples=1 evals=1
println("Undecided domain: ", round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")

# Save 

filename = splitext(PROGRAM_FILE)[1]
if allow_exists_and_forall_bisection
    outfile = "$(filename)_11_$(ϵ_x)_$(ϵ_p)_refined.png"
elseif allow_exists_or_forall_bisection
    outfile = "$(filename)_11_$(ϵ_x)_$(ϵ_p)_subdivided.png"
else
    outfile = "$(filename)_11_$(ϵ_x).png"
end

save_drawing(X_0, inn, out, delta, outfile)