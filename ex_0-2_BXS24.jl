using NeuralVerification
using IntervalArithmetic

include("utils.jl")
include("pave.jl")

nnet = read_nnet("BXS24.nnet"; last_layer_activation = NeuralVerification.Id())

n = 1
p = 3
f_fun = [x -> NeuralVerification.compute_output(nnet, [x[1:2]...])[1] - NeuralVerification.compute_output(nnet, [x[1:2]...])[2] - x[3]]
Df_fun = [x -> [(get_gradient(nnet, x[1:2])[1, :] - get_gradient(nnet, x[1:2])[2, :])... -1]]
problem = Problem(f_fun, Df_fun)
qvs = []
qcp = QuantifiedConstraintProblem(problem, qvs, [qvs], p, n)
X_0 = IntervalBox(interval(0, 20), interval(0, 20))
p_in = []
p_out = deepcopy(p_in)
G = [interval(0, 100)]

ϵ_x = 1.0
ϵ_p = nothing
inn, out, delta = pave_11(X_0, p_in, p_out, G, qcp, ϵ_x, ϵ_p, false, false)

p = plot()
draw(p, X_0, inn, out, delta)

savefig("0-2_BXS24.png")
