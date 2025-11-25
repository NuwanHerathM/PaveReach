# Detailed explanation

We explain here the core of `ex_5-1_running_example.jl`.

The example in `ex_5-1_running_example.jl` sets up the scalar problem:
$$\{ x \in [-5, 5] \, | \, \forall p_1 \in [0, 1/4], f(x, p_1) = p_1^2 - (x - 1)(x - 2)(x - 3) \in [-1/4, 1/4] \}$$

- Domain for $x$: $X_0 = [-5, 5]$
- Parameters and quantifiers: $∀ p1 ∈ [0, 1/4]$
- Target: $\mathbb{G} = [-1/4, 1/4]$

Command line arguments allow the user to choose: 
- $\mathcal{O}^{IN}$,
- $\mathcal{O}^{OUT}$
- $\epsilon_\mathbb{X}$
- and $\epsilon_\mathbb{P}$.

Additionally, one can choose the method to subdivide the parameters:
- `allow_exists_or_forall_bisection` for Section 4.1,
- `allow_exists_and_forall_bisection` for Section 4.2.
They are set with command line options.

The script:
1. Defines the problem.
    ```julia
    n = 1
    p = 3
    @variables x[1:p]
    f_num = [x[2]^2-(x[1]-1)*(x[1]-2)*(x[1]-3)-x[3]]
    f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
    problem = Problem(f_fun, Df_fun) 
    qcp = QuantifiedConstraintProblem(problem, [(Forall, 2)], [[(Forall, 2)]], p, n)
    ```
2. $\mathbb{X}_0$, $\mathbb{P}$ and $\mathbb{G}$ are treated differently.
    * $\mathbb{X}_0$ is an `IntervalBox`.
    ```julia
    X_0 = IntervalBox(interval(0, 5))
    ```
    * $\mathbb{P}$ is a list of list with one `Interval` for each component, by anticipation of the bisections to come: if parameters are subdivided, `[[p_1], [p_2], ...]` will become `[[p_1^1, p_1^2, ...], [p_2^1, ...], ...]`.
    ```julia
    p_in = [[interval(0, 1/4)]]
    p_out = deepcopy(p_in)
    ```
    * $\mathbb{G}$ is just a list of intervals, since it will not be bisected: `[g_1, g_2, ...]`.
    ```julia
    G = [interval(-1/4, 1/4)]
    ```
3. If the method of Section 4.2 has been chosen, modifies the parameters accordingly.
    ```julia
    if allow_exists_and_forall_bisection
        refine_in!(p_in, qcp.qvs, ϵ_p, qcp.p, qcp.n)
        refine_out!(p_out, qcp.qvs, ϵ_p, qcp.p, qcp.n)
    end
    ```
4. Runs the chosen paving function. Times it using `@btime` from the `BenchmarkTools` package.
    ```julia
    @btime (global inn, out, delta = paving_function(X_0, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection))
    ```
5. Prints the proportion of undecided domain.
    ```julia
    println("Undecided domain: ", round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")
    ```