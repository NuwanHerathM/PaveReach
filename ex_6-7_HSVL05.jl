include("pave2.jl")

n = 1
p = 4
@variables x[1:p]
f_num = [x[1]^2+x[2]^2+2*x[1]*x[2]-20*x[1]-20x[2]+100-x[3]-x[4]]
f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
problem = Problem(f_fun, Df_fun) 
qvs = [(Forall, 3), (Exists, 2)]
qe = QuantifierProblem(problem, qvs, [qvs], p, n)
X_0 = IntervalBox(interval(0, 6))
p_in = [[interval(2, 8)], [interval(6, 8)]]
p_out = deepcopy(p_in)
G = [[interval(0, 0)]]

is_refined = false
println(if is_refined "Refined" else "Not refined" end)
ϵ = 0.01
println("ϵ  = ", ϵ)
allow_normal_p_bisect = true
println(if allow_normal_p_bisect "Normal bisection on P" else "No standard bisection on P" end)
if is_refined
    refine_in!(p_in, qe.qvs, ϵ, qe.p, qe.n)
    refine_out!(p_out, qe.qvs, ϵ, qe.p, qe.n)
end
println(length.(p_in))
println(length.(p_out))
inn, out, delta = pave(X_0, p_in, p_out, G, qe, ϵ, is_refined, allow_normal_p_bisect)

println("Number of elements of inn: ", merge_intervals(first.(inn)))
println("Number of elements of out: ", merge_intervals(first.(out)))
println("Number of elements of delta: ", merge_intervals(first.(delta)))

println(round(volume_boxes(delta)/volume_box(X_0)*100, digits=1), " %")

draw(X_0, [inn], [out], [delta])
gui()