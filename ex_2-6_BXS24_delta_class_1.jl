include("utils.jl")
include("pave.jl")

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
allow_exists_and_forall_bisection = parsed_args["refine"]
allow_exists_or_forall_bisection = parsed_args["subdivide"]

ϵ_x = [0.1]
ϵ_p = [1.0, 1.0, 1.0, 1.0]

if isnothing(ϵ_p)
    @assert !allow_exists_and_forall_bisection "eps_p was not provided. Refinement requires eps_p. Use --help for more information."
    @assert !allow_exists_or_forall_bisection "eps_p was not provided. Subdivision requires eps_p. Use --help for more information."
else
    if !allow_exists_and_forall_bisection && !allow_exists_or_forall_bisection
        @warn "eps_p was provided, but will not be used. Provide options --refine or --subdivide in order to use eps_p."
    end
end
@assert nand(allow_exists_and_forall_bisection, allow_exists_or_forall_bisection) "Refinement and subdivision are mutually exclusive. Use --help for more information."
# ------------------------------------------------------

filename = splitext(PROGRAM_FILE)[1]

nnet = read_nnet("BXS24.nnet"; last_layer_activation = NeuralVerification.Id())
N(x) = NeuralVerification.compute_output(nnet, x)
DN(x) = get_gradient(nnet, x)

confidence(vect::AbstractVector) = vect[1] - vect[2]
confidence(mat::Matrix) = mat[1, :] - mat[2, :]

function perturbation(x)
    return [x[1] + x[3], x[2] + x[4]]
end

function gradient_perturbation(x)
    return [1 0 1 0; 0 1 0 1]
end

n = 2
p = 7
f_fun = [x -> (confidence(N(x[2:3])) - x[1] - x[6]), x -> (confidence(N(perturbation(x[2:5]))) - x[7])]
Df_fun = [x -> [-1 confidence(DN(x[2:3]))... 0 0 -1 0],
        x -> [0 confidence(DN(perturbation(x[2:5])) * gradient_perturbation(x[2:5]))... 0 -1]]
problem = Problem(f_fun, Df_fun, [[1], [2]])
qvs = [(Forall, 2), (Forall, 3), (Forall, 4), (Forall, 5)]
qcp = QuantifiedConstraintProblem(problem, qvs, [qvs, qvs], p, n)
X_0 = IntervalBox(interval(0, 10))
ϵ_max = 0.25
p_in = [[interval(-10, 10)], [interval(-10,-10)], [interval(-ϵ_max, ϵ_max)], [interval(-ϵ_max, ϵ_max)]]
p_out = deepcopy(p_in)
G = [interval(minus_inf, -0.0001), interval(0.0001, plus_inf)]

if allow_exists_and_forall_bisection
    refine_in!(p_in, qcp.qvs, ϵ_p, qcp.p, qcp.n)
    refine_out!(p_out, qcp.qvs, ϵ_p, qcp.p, qcp.n)
end

inn, out, delta = pave_12(X_0, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection)

if parsed_args["save"]
    p = plot()
    draw(p, X_0, inn, out, delta)

    plot(p, size = (600, 80), grid = false)

    if allow_exists_and_forall_bisection
        outfile = "$(filename)_$(ϵ_x)_$(ϵ_p)_refined.png"
    elseif allow_exists_or_forall_bisection
        outfile = "$(filename)_$(ϵ_x)_$(ϵ_p)_subdivided.png"
    else
        outfile = "$(filename)_$(ϵ_x).png"
    end
    savefig(outfile)
    println("The result was saved in $(outfile).")
end
