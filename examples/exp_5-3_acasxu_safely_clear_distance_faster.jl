include("../src/PaveReach.jl")

using BenchmarkTools

# Load the neural network

# ACAS Xu model from the Reluplex paper
# It is from 2017
# https://github.com/guykatzz/ReluplexCav2017/tree/master/nnet
modelpath = joinpath(@__DIR__, "models/ACASXU_run2a_1_1_batch_2000.nnet")
nnet = read_nnet(modelpath)
N(x) = NeuralVerification.compute_output(nnet, x)
DN(x) = get_gradient(nnet, x)
dN(i, x) = get_partial_derivative(nnet, i, x)
N(x::Tuple) = N(collect(x))
DN(x::Tuple) = DN(collect(x))
dN(i, x::Tuple) = dN(i, collect(x))

const mins = [0.0,-3.141593,-3.141593,100.0,0.0]
const maxes = [60760.0,3.141593,3.141593,1200.0,1200.0]
const means = [19791.091,0.0,0.0,650.0,600.0]
const ranges = [60261.0,6.28318530718,6.28318530718,1100.0,1200.0]
const output_mean = 7.5188840201005975
const output_range = 373.94992

function normalize_input(index, value)
    min = mins[index]
    max = maxes[index]
    mean = means[index]
    range = ranges[index]

    if value < min
        value = min
    elseif value > max
        value = max
    end

    return (value - mean) / range
end

function normalize_input(index, value::IntervalArithmetic.Interval{T}) where T <: AbstractFloat
    min = mins[index]
    max = maxes[index]
    mean = means[index]
    range = ranges[index]

    lo = 
    if value.hi < min
        value = interval(min)
    elseif value.lo > max
        value = interval(max)
    end

    return (value - mean) / range
end

function normalize_output(value)
    return (value - output_mean) / output_range
end

# Set the problem
# ρ ∈ [0, 10000] ∣ N_1(ρ, 0, 0, 1000, 600) ≤ 22 ⟹ N_1(ρ, 0, 0, 1000, 600 - 5) ≤ 22.5

p = 1
θ = 0.0
ψ = 0.0
v_o = 1000.0
v_i = 600.0
x_constant_1 = [normalize_input(2, interval(θ)), normalize_input(3, interval(ψ)), normalize_input(4, interval(v_o)), normalize_input(5, interval(v_i))]
x_constant_2 = [normalize_input(2, interval(θ)), normalize_input(3, interval(ψ)), normalize_input(4, interval(v_o)), normalize_input(5, interval(v_i - 5))]
f_fun = Function[(x...) -> N([x[1], x_constant_1...])[1], (x...) -> N([x[1], x_constant_2...])[1]]
Df_fun = ColJacobianMatrix([(x...) -> [dN(1, [x[1], x_constant_1...])[1], dN(1, [x[1], x_constant_2...])[1]]])
sizes = [1, 1]
formula = parseformula("(¬ P_1) ∨ P_2")
qvs = QuantifiedVariable[]
qvs_relaxed = [qvs]
parameters = ProblemParameters(formula, f_fun, Df_fun, sizes, qvs, qvs_relaxed, p)
X = interval(0.0, 10000.0)
X_0 = normalize_input(1, X)
P = []
δ_1 = 22
δ_2 = 22.5
G = [interval(minus_inf, normalize_output(δ_1)), interval(minus_inf, normalize_output(δ_2))]
domains = ProblemDomains(P, G)

ϵ_x = [0.0001]
configuration = PavingConfiguration(ϵ_x)

# Pave

verbose = false
benchmark = @benchmarkable (global inn, out, delta = pave(X_0, parameters, domains, configuration, verbose)) samples=1 evals=1
results = run(benchmark)
timing = round(median(results).time/1e6, digits=3) # time/1e6 -> timing in milliseconds

# Save the timing and the volumes in a log file

filename = splitext(basename(PROGRAM_FILE))[1]
base = "$(filename)_$(ϵ_x)"

logfile = base * ".log"
save_log(timing, inn, out, delta, ranges[1], @__DIR__, logfile)