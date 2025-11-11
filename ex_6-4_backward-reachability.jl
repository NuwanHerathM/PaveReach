include("pave2.jl")

filename = splitext(PROGRAM_FILE)[1]

n = 1
p = 9
@variables x[1:p]
f_num = [(-0.5+0.1*x[1]+(1+0.01*x[4])*x[8]+1.31*10^(-7)*x[3]*x[8]^2)^2 + (-0.05+0.1*x[2]+(0.01*x[6]+0.01*x[7]*x[8])*x[8]+0.005*x[5]*x[8]^2)^2 - x[9]]
f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
problem = Problem(f_fun, Df_fun)
qvs = [(Exists, 7), (Exists, 6), (Forall, 4), (Exists, 3), (Exists, 5), (Exists, 8)]
qe = QuantifierProblem(problem, qvs, [qvs], p, n)
X_0 = IntervalBox(interval(-2, 12), interval(-7, 7))
p_in = [[interval(-1, 1)], [interval(-1, 1)], [interval(-1, 1)], [interval(-1, 1)], [interval(-1, 1)], [interval(-1, 1)], [interval(0, 0.5)]]
p_out = deepcopy(p_in)
G = [interval(0, 1)]

is_refined = false
println(if is_refined "Refined" else "Not refined" end)
ϵ = 0.1
println("ϵ  = ", ϵ)
allow_normal_p_bisect = false
println(if allow_normal_p_bisect "Normal bisection on P" else "No standard bisection on P" end)
if is_refined
    refine_in!(p_in, qe.qvs, ϵ, qe.p, qe.n)
    refine_out!(p_out, qe.qvs, ϵ, qe.p, qe.n)
end

@assert nand(is_refined, allow_normal_p_bisect)

factor = 5.0

inn, out, delta = pave_11(X_0, p_in, p_out, G, qe, ϵ, factor, is_refined, allow_normal_p_bisect)
p11 = plot()
draw(p11, X_0, inn, out, delta)
println(round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")

str_refined = is_refined ? "_refined" : ""
str_subdivided = allow_normal_p_bisect ? "_subdivided" : ""

savefig("$(filename)_$(ϵ)$(str_refined)$(str_subdivided).png")
gui()