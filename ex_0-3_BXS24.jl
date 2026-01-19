using NeuralVerification
using IntervalArithmetic

include("utils.jl")
include("pave.jl")

nnet = read_nnet("BXS24.nnet"; last_layer_activation = NeuralVerification.Id())
N(x) = NeuralVerification.compute_output(nnet, x)
DN(x) = get_gradient(nnet, x)

confidence(vect::AbstractVector) = vect[1] - vect[2]
confidence(mat::Matrix) = mat[1, :] - mat[2, :]

function perturbation(x)
    return [x[1] + x[3], x[2] + x[4]]
end

function gradient_perturbation(x)
    return [x[1] 0 x[3] 0; 0 x[2] 0 x[4]]
end

δ = 0

n = 1
p = 5
f_fun = [x -> (confidence(N(perturbation(x[1:4]))) - x[5])]
ϕ = [x -> confidence(N(perturbation(x[1:4])))]
Df_fun = [x -> [confidence(DN(perturbation(x[1:4])) * gradient_perturbation(x[1:4]))... -1]]
problem = Problem(f_fun, Df_fun, ϕ, [[1]])
qvs = [(Forall, 3), (Forall, 4)]
qcp = QuantifiedConstraintProblem(problem, qvs, [qvs], p, n)
X_0 = IntervalBox(interval(0, 20), interval(0, 20))
p_in = [interval(-0.5, 0.5), interval(-0.5, 0.5)]
p_out = deepcopy(p_in)
G = [interval(δ, 100)]

ϵ_x = 0.5
ϵ_p = nothing
inn, out, delta = pave_12(X_0, p_in, p_out, G, qcp, ϵ_x, ϵ_p, false, false)

p = plot()
draw(p, X_0, inn, out, delta)

savefig("0-3_BXS24_$(δ)_$(ϵ_x).png")
