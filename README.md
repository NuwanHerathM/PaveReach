# PaveReach

## Quick start

Try the running example from the article (Example 3.1):
```julia
# Non-interactive run (prints inn/out/delta)
julia 3-2_runningexample.jl

# Run and keep REPL open for interactive plotting / display
julia -i 3-2_runningexample.jl -d
```
Try a 2D toy example:
```julia
# Non-interactive run (prints size of inn/out/delta)
julia 0-1_disk.jl

# Run and keep REPL open for interactive plotting / display
julia -i 0-1_disk.jl -d
```

## What the examples do

### Detailed explanation for a 1D example (`3-2_runningexample.jl`)

The example in `3-2_runningexample.jl` sets up the scalar problem:
$$\{ x \in [-5, 5] \, | \, \forall p_1 \in [0, 1/4], \exists z \in [-1/4, 1/4], f(x, p_1, z) = p_1^2 - (x - 1)(x - 2)(x - 3) - z = 0 \}$$

- Quantifiers: $∀ p1 ∈ [0, 1/4], ∃ z ∈ [-1/4, 1/4]$
- Domain for x: $X_0 = [-5, 5]$
- Paving precision: $\epsilon = 0.1$

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
# problem := f, Df
problem = Problem(f_fun, Df_fun)
# problem, [∀ p1, ∃ z], p, n
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

### Brief explanation for a 2D example (`0-1_disk.jl`)

The example in `0-1_disk.jl` sets up the scalar problem:
$$\{ x \in [-5, 5] \times [-5, 5] \, | \, \exists z \in [0, 16], f(x, z) = x_1^2 + x_2^2 - z = 0 \}$$

- Quantifier: $∃ z ∈ [0, 16]$
- Domain for x: $X_0 = [-5, 5] \times [-5, 5]$
- Paving precision: $\epsilon = 0.1$

The script:
1. Constructs a `QuantifierProblem`. Beware that the free variable $x$ is affected to the variables `x[1]` and `x[2]`. By default, $z$ is then affected to `x[3]`.
```julia
# 1
n = 1
# Dimension of the domain of f
p = 3
# x[1] := x_1, x[2] := x_2, x[3] := z
@variables x[1:p]
# f(x, p_1, z) = x_1^2 + x_2^2 - z
f_num = [x[1]^2 + x[2]^2 - x[3]]
f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
# problem := f, Df
problem = Problem(f_fun, Df_fun)
# f, Df, [∃ z], p, n
qe = QuantifierProblem(problem, [(Exists, 3)], p, n)
```
2. Builds initial inner/outer refinement trees and runs `pave_11`.
```julia
# Domain to be paved X_0
X_0 = IntervalBox(interval(-5, 5), interval(-5, 5))
# z ∈ Z = [0, 16]
intervals = [interval(0, 16)]
# Precision
eps = 0.1
# Paving using the 1st criterion for IN and the 1st criterion for OUT
pz_in_0, pz_out_0 = make_pz_11(intervals, qe)
inn, out, delta = pave_11(pz_in_0, pz_out_0, qe, X_0, eps)
```
3. Same as the previous example...
4. ...

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

## Parameter refinement (P, Z)

No refinement:
```julia
inn, out, delta = pave_11(pz_in_0, pz_out_0, qe, X_0, eps, is_refined=false)
```
Refinement:
```julia
inn, out, delta = pave_11(pz_in_0, pz_out_0, qe, X_0, eps)
```
or
```julia
inn, out, delta = pave_11(pz_in_0, pz_out_0, qe, X_0, eps, is_refined=true)
```

## Conjunction and disjunction

For a conjunction
$$\{ \dots \, | \, \dots, f_1(x, p, z) = 0  \, \wedge \, f_2(x, p, z) = 0 \, \wedge \, f_3(x, p, z) = 0 \}$$
```julia
f_fun_1, Df_fun_1 = ...
problem_1 = Problem(f_fun_1, Df_fun_1)
f_fun_2, Df_fun_2 = ...
problem_2 = Problem(f_fun_2, Df_fun_2)
f_fun_3, Df_fun_3 = ...
problem_3 = Problem(f_fun_3, Df_fun_3)
# Conjunction
problem = AndProblem([problem_1, problem_2, problem_3])
```

For a combination of disjunction and conjunction
$$\{ \dots \, | \, \dots, (f_1(x, p, z) = 0  \, \vee \, f_2(x, p, z) = 0) \, \wedge \, f_3(x, p, z) = 0 \}$$
```julia
f_fun_1, Df_fun_1 = ...
problem_1 = Problem(f_fun_1, Df_fun_1)
f_fun_2, Df_fun_2 = ...
problem_2 = Problem(f_fun_2, Df_fun_2)
f_fun_3, Df_fun_3 = ...
problem_3 = Problem(f_fun_3, Df_fun_3)
# Disjunction -- subproblem
sub_problem = OrProblem([problem_1, problem_2])
# Conjunction -- problem
problem = AndProblem([sub_problem, problem_3])
```

See `0-3_disk_intersection.jl` for an example of a conjuction.