include("../src/PaveReach.jl")

using BenchmarkTools

# Load the neural network

# ACAS Xu model from the Reluplex paper
# It is from 2017
# https://github.com/guykatzz/ReluplexCav2017/tree/master/nnet
nnet = read_nnet("models/ACASXU_run2a_1_1_batch_2000.nnet")
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

# ϕ_1 coc < 1500

# Pave ρ, θ, v_o and v_i
p = 4
x_constant = [normalize_input(3, interval(0.0))]
f_fun = Function[(x...) -> N([x[1], x[2], x_constant..., x[3], x[4]])[1]]
Df_fun = ColJacobianMatrix([(x...) -> [dN(1, [x[1], x[2], x_constant..., x[3], x[4]])[1]],
            (x...) -> [dN(2, [x[1], x[2], x_constant..., x[3], x[4]])[1]],
            (x...) -> [dN(4, [x[1], x[2], x_constant..., x[3], x[4]])[1]],
            (x...) -> [dN(5, [x[1], x[2], x_constant..., x[3], x[4]])[1]]])
sizes = [1]
formula = parseformula("P_1")
qvs = QuantifiedVariable[]
qvs_relaxed = [qvs]
parameters = ProblemParameters(formula, f_fun, Df_fun, sizes, qvs, qvs_relaxed, p)
X = [interval(55947.0, 60760.0), interval(-π, π), interval(1145.0, 1200.0), interval(0.0, 60.0)]
X_0 = [normalize_input(1, X[1]), normalize_input(2, X[2]), normalize_input(4, X[3]), normalize_input(5, X[4])]
P = []
G = [interval(minus_inf, normalize_output(1500.0) - strict_epsilon)]
domains = ProblemDomains(P, G)

ϵ_x = [0.01, 0.005, 0.01, 0.01]
configuration = PavingConfiguration(ϵ_x)

# Pave

verbose = false
benchmark = @benchmarkable (global inn, out, delta = pave(X_0, parameters, domains, configuration, verbose)) samples=1 evals=1
results = run(benchmark)
timing = round(median(results).time/1e6, digits=3) # time/1e6 -> timing in milliseconds

# Save paving information in a .log file

filename = splitext(PROGRAM_FILE)[1]
base = "$(filename)_$(ϵ_x)"
logfile = base * ".log"

save_log(timing, inn, out, delta, logfile)
