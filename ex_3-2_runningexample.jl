include("pave2.jl")

n = 1
p = 3
@variables x[1:p]
f_num = [x[2]^2-(x[1]-1)*(x[1]-2)*(x[1]-3)-x[3]]
f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
problem = Problem(f_fun, Df_fun) 
qe = QuantifierProblem(problem, [(Forall, 2)], [[(Forall, 2)]], p, n)
X_0 = IntervalBox(interval(0, 5))
p_in = [[interval(0, 1/4)]]
p_out = deepcopy(p_in)
G = [[interval(-1/4, 1/4)]]

is_refined = false
println(if is_refined "Refined" else "Not refined" end)
ϵ = 0.1
println("ϵ  = ", ϵ)
allow_normal_p_bisect = true
println(if allow_normal_p_bisect "Normal bisection on P" else "No standard bisection on P" end)
if is_refined
    refine_in!(p_in, qe.qvs, ϵ, qe.p, qe.n)
    refine_out!(p_out, qe.qvs, ϵ, qe.p, qe.n)
end
inn, out, delta = pave(X_0, p_in, p_out, G, qe, ϵ, is_refined, allow_normal_p_bisect)

println("Number of elements of inn: ", length(inn))
println("Number of elements of out: ", length(out))
println("Number of elements of delta: ", length(delta))

println(round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")

draw(X_0, [inn], [out], [delta])
gui()