include("utils.jl")
include("pave.jl")
include("read_onnx.jl")

# ------------------------------------------------------
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "delta"
            help = "confidence delta"
            arg_type = Float64
            required = true
        "eps_x"
            help = "epsilon for the paving"
            arg_type = Float64
            required = true
        "eps_p"
            help = "epsilon for the parameters"
            arg_type = Float64
            required = false
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
δ = parsed_args["delta"]

ϵ_x = parsed_args["eps_x"]
ϵ_p = parsed_args["eps_p"]
allow_exists_and_forall_bisection = parsed_args["refine"]
allow_exists_or_forall_bisection = parsed_args["subdivide"]


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

nnet = read_onnx_mlp("normed_etcs_2x250_nodes.onnx")
N(x) = NeuralVerification.compute_output(nnet, x)
DN(x) = get_gradient(nnet, x)

confidence(vect::AbstractVector) = vect[1] - vect[2]
confidence(mat::Matrix) = mat[1, :] - mat[2, :]

gradient_variables = [1 0; 0 1; 0 0]

x_r = 1

n = 1
p = 3
f_fun = [x -> (confidence(N([x[1], x[2], x_r])) - x[3])]
Df_fun = [x -> [confidence(DN(IntervalBox([x[1], x[2], x_r])) * gradient_variables)... -1]]
problem = Problem(f_fun, Df_fun, [[1]])
qvs = []
qcp = QuantifiedConstraintProblem(problem, qvs, [qvs], p, n)
X_0 = IntervalBox(interval(0, 1), interval(0, 1))
p_in = []
p_out = deepcopy(p_in)
G = [interval(δ, plus_inf)]

inn, out, delta = pave_12(X_0, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection)

if parsed_args["save"]
    p = plot()
    draw(p, X_0, inn, out, delta)

    plot(p)

    if allow_exists_and_forall_bisection
        outfile = "$(filename)_$(δ)_$(ϵ_x)_$(ϵ_p)_refined.png"
    elseif allow_exists_or_forall_bisection
        outfile = "$(filename)_$(δ)_$(ϵ_x)_$(ϵ_p)_subdivided.png"
    else
        outfile = "$(filename)_$(δ)_$(ϵ_x).png"
    end
    savefig(outfile)
    println("The result was saved in $(outfile).")
end
