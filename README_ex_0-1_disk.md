# Detailed explanation for a 2D example (`ex_0-1_disk.jl`)

> Note: the example is degenerated, there is no parameter $\mathbb{P}$.

The example in `0-1_disk.jl` sets up the scalar problem:
$$\{ x \in [-5, 5] \times [-5, 5] \, | \, \exists z \in [0, 16], f(x, z) = x_1^2 + x_2^2 - z = 0 \}.$$

- Domain for x: $\mathbb{X}_0 = [-5, 5] \times [-5, 5]$
- Parameters and quantifiers: none
- Target: $\mathbb{G} = [0, 16]$
- Paving precision: $\epsilon_{\mathbb{X}} = 0.1$
- Parameter subdivision precision: $\epsilon_{\mathbb{P}} = 0.5$ (will not be used since there is no $\mathbb{P}$)

The script:
1. Constructs a `QuantifiedConstraintProblem`.
    * Defines the dimensions of the domain and the codomain of $f: \mathbb{R}^p \to \mathbb{R}^n$.
        ```julia
        # Dimension of the codomain of f
        n = 1
        # Dimension of the domain of f
        p = 3
        ```
    * Writes the symbolic expression of $f: (x,p) \mapsto f(x,p)$. Beware that $x$ is affected to the first variables, `x[1]` and `x[2]`, and that $z$ is then affected to the last $n$ variables, here only `x[3]`. If there was a $p$ it would be affected to the remaining variables.
        ```julia
        # x[1] := x_1, x[2] := x_2, x[3] := z
        @variables x[1:p]
        # f(x, p_1, z) = x_1^2 + x_2^2 - z
        f_num = [x[1]^2 + x[2]^2 - x[3]]
        ```
    * Builds the function $f$ and its Jacobian.
        ```julia
        # Symbolic expression to functions
        f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
        ```
    * Creates the quantified problem.
        ```julia
        # problem := f, Df
        problem = Problem(f_fun, Df_fun)
		# Quantified variables
        qvs = []
        # Quantifier elimination problem
		qcp = QuantifiedConstraintProblem(problem, qvs, [qvs], p, n)
        ```
    * Defines the domain of the parameters $\mathbb{P}$.
        ```julia
        # Parameters
        p_in = []
        p_out = deepcopy(p_in)
        ```
    * Defines the target domain $\mathbb{G} \in \mathbb{R}^n$.
        ```julia
        # Target
        G = [interval(0, 16)]
        ```
2. Sets the precision and chooses a method to subdivide the parameters.
	```julia
	# Precision
	ϵ_x = 0.1
	ϵ_p = 0.5
	# Subdivision as in Section 4.1
	allow_exists_or_forall_bisection = true
	# Subdivision as in Section 4.2
	allow_exists_and_forall_bisection = false
	```
	Checks that both methods are not set to `true`.
	```julia
	# Assert both subdivision methods have not been chosen
	@assert nand(allow_exists_and_forall_bisection, allow_exists_or_forall_bisection)
	```
	If the method of Section 4.2 has been chosen, modifies the parameters accordingly.
	```julia
	# Refine (points/subdivision)
	if allow_exists_and_forall_bisection
		refine_in!(p_in, qcp.qvs, ϵ_p, qcp.p, qcp.n)
		refine_out!(p_out, qcp.qvs, ϵ_p, qcp.p, qcp.n)
	end
	```
3. Prints a recap of the choices.
	```julia
	# Print info
	println("ϵ_x  = ", ϵ_x)
	println("ϵ_p  = ", ϵ_p)
	println(if allow_exists_and_forall_bisection "Refined" else "Not refined" end)
	println(if allow_exists_or_forall_bisection "Normal bisection on P" else "No standard bisection on P" end)
	```
4. Runs `pave_11`.
	```julia
	# Paving with O^IN using P and G, as well as O^OUT using P and G
	inn, out, delta = pave_11(X_0, p_in, p_out, G, qcp, ϵ_x, ϵ_p, allow_exists_and_forall_bisection, allow_exists_or_forall_bisection)
	```
5. Prints the proportion of undecided domain.
	```julia
	# Percentage of undecided domain
	println("Undecided domain: ", round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")
	```
6. Creates a plot and saves it in an well-named file.
	```julia
	# Construct the plot
	p = plot()
	draw(p, X_0, inn, out, delta)
	plot(p)
	# Name the output file according the method
	if allow_exists_and_forall_bisection
		outfile = "$(filename)_11_$(ϵ_x)_$(ϵ_p)_refined.png"
	elseif allow_exists_or_forall_bisection
		outfile = "$(filename)_11_$(ϵ_x)_$(ϵ_p)_subdivided.png"
	else
		outfile = "$(filename)_11_$(ϵ_x).png"
	end
	# Save the output file
	savefig(outfile)
	println("The result was saved in $(outfile).")
	```