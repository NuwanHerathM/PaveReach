# PaveReach

## Quick start

Try the running example from the article (Example 3.1):
```julia
# Non-interactive run (prints inn/out/delta)
julia ex_runningexample.jl

# Run and keep REPL open for interactive plotting / display
julia -i ex_runningexample.jl -d
```

## What the example does

The example in `ex_runningexample.jl` sets up the scalar problem:
$$\{ x \in [-5, 5] \, | \, \forall p_1 \in [0, 1/4], \exists z \in [-1/4, 1/4], f(x, p_1, z) = p_1^2 - (x - 1)(x - 2)(x - 3) - z\}$$

- Quantifiers: ∀ p1 ∈ [0, 1/4], ∃ z ∈ [-1/4, 1/4]
- Domain for x: X_0 = [-5, 5]
- Paving precision: `eps = 0.1`

The script:
1. Constructs a `QuantifierProblem`. Beware that the free variable $x$ is affected to the variable `x[1]`. The order does not matter for $p_1$ and $z$ however.
```julia
# 1
n = 1
# Dimension of the domain of f
p = 3
# x[1] := x, x[2] := p_1, x[3] := z
@variables x[1:p]
# f(x, p_1, z) = p_1^2 - (x - 1)(x - 2)(x - 3) - z
f_num = [x[2]^2-(x[1]-1)*(x[1]-2)*(x[1]-3)-x[3]]
f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
# f, Df, [∀ p1, ∃ z], p, n
qe = QuantifierProblem(f_fun, Df_fun, [(Forall, 2), (Exists, 3)], p, n)
```
2. Builds initial inner/outer refinement trees and runs `pave_11`.
```julia
# Domain to be paved X_0
X_0 = IntervalBox(interval(-5, 5))
# p_1 ∈ [0, 1/4], z ∈ Z
intervals = [interval(0, 1/4), interval(-1/4, 1/4)]                                               
# Precision
eps = 0.1
# Paving using the 1st criterion for IN and the 1st criterion for OUT
pz_in_0, pz_out_0 = make_pz_11(intervals, qe)
inn, out, delta = pave_11(pz_in_0, pz_out_0, qe, X_0, eps)
```
3. Prints the results: `inn`, `out`, and `delta` (inner approximation, outer approximation, and undecided / boundary set).
```julia
print_inn_out_delta(inn, out, delta)
```
4. Optionally displays a visualization when the `--display` (`-d`) flag is given.
After running, the script prints three objects:
- `inn` — green
- `out` — blue
- `delta` — yellow

## Oracles

To choose the 1st criterion for IN and the 1st criterion for OUT, select `make_pz_11` and `pave_11`
```julia
pz_in_0, pz_out_0 = make_pz_11(intervals, qe)
inn, out, delta = pave_11(pz_in_0, pz_out_0, qe, X_0, eps)
```
To choose the 1st criterion for IN and the 2st criterion for OUT, select `make_pz_12` and `pave_12`
```julia
pz_in_0, pz_out_0 = make_pz_12(intervals, qe)
inn, out, delta = pave_12(pz_in_0, pz_out_0, qe, X_0, eps)
```
And so on...