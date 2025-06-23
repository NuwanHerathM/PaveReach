#!/usr/local/bin/julia

using BenchmarkTools
using IntervalArithmetic

include("genreach2.jl")
include("cell.jl")

@variables x[1:3]
QE5=([x[1]^2/4.0+(x[2]+1.0)*(x[3]+2.0)+(x[3]+3.0)^2],["exists", 1, "forall", 2, "exists", 3],3,1)
# print_QE(QE5)

# println("Results for ex5 -order 0")
# res_QE5_o0=QEapprox_o0(QE5[1], QE5[2], [["exists", 1, "forall", 2, "exists", 3]], QE5[3], QE5[4], [interval(-1,1), interval(-1,1), interval(-1,1)])
# println(res_QE5_o0)
# res_QE5_o0=QEapprox_o0(QE5[1], QE5[2], [["exists", 1, "forall", 2, "exists", 3]], QE5[3], QE5[4], [interval(-1,1), interval(-1,0), interval(-1,1)])
# println(res_QE5_o0)
# res_QE5_o0=QEapprox_o0(QE5[1], QE5[2], [["exists", 1, "forall", 2, "exists", 3]], QE5[3], QE5[4], [interval(-1,1), interval(0,1), interval(-1,1)])
# println(res_QE5_o0)

function approximate(box::IntervalBox)
    return QEapprox_o0(QE5[1], QE5[2], [["exists", 1, "forall", 2, "exists", 3]], QE5[3], QE5[4], [box[i] for i=1:length(box)])
end

function get_quantifier_from_QE(QE::Tuple{Vector{Num}, Vector{Any}, Int64, Int64}, dim::Int)
    i = 2
    while QE[2][i] != dim
        i += 2
    end
    return QE[2][i-1]
end

function update!(cell::ApproximationCell)
    cell_1 = cell.children[1]
    cell_2 = cell.children[2]
    dim = 1
    while cell_1.box[dim] == cell_2.box[dim]
        dim += 1
    end
    quantifier = get_quantifier_from_QE(QE5, dim)
    if quantifier == "exists"
        cell.inner = cell_1.inner ∪ cell_2.inner
        cell.outer = cell_1.outer ∪ cell_2.outer
    else
        cell.inner = cell_1.inner ∩ cell_2.inner
        cell.outer = cell_1.outer ∩ cell_2.outer
    end

    if !isroot(cell)
        update!(cell.parent)
    end
end

################################

eps_tolerance = 0.25
eps_diam = 0.2

@time begin

box = IntervalBox(interval(-1,1), interval(-1,1), interval(-1,1)) 

inner, outer = approximate(box)

cell_0 = make_approximationcell_root(box, inner[1].dat, outer[1].dat)
root = cell_0

# println("Initial tolerance: ", tolerance(root))
# println(root)


nb_bisections = 0
l = [root]
while tolerance(root) > eps_tolerance

    cell = popfirst!(l)

    if diam(cell.box) <= eps_diam
        break
    end

    box_1, box_2 = bisect(cell.box, 0.5)
    global nb_bisections += 1

    inner_1, outer_1 = approximate(box_1)
    cell_1 = make_approximationcell_leaf(box_1, inner_1[1].dat, outer_1[1].dat, cell)
    
    inner_2, outer_2 = approximate(box_2)
    cell_2 = make_approximationcell_leaf(box_2, inner_2[1].dat, outer_2[1].dat, cell)

    push!(l, cell_1)
    push!(l, cell_2)

    update!(cell)

end

end

println("Tolerance: ", tolerance(root))
println(root)
println("Number of bisections: ", nb_bisections)

################################

# using Plots; pythonplot()
# xlims!((-1,1))
# ylims!((-1,1))

# rectangle(p, q) = Shape([p[1],q[1],q[1],p[1]], [p[2],p[2],q[2],q[2]])

# current_plot = plot()

# function draw(cell::Cell)
#     if isleaf(cell)
#         box = cell.box
#         xs = (box[1].lo, box[2].lo)
#         ys = (box[1].hi, box[2].hi)
#         plot!(rectangle(xs, ys), color=:blue)
#     else
#         for child in cell.children
#             draw(child)
#         end
#     end
# end

# draw(cell_0)
# plot(current_plot, legend=:false)
# gui()

# dump(stripeset.stripes[2])

# cell_1 = Cell(box_1)

# g = x[1]^2/4.0+(x[2]+1.0)*(x[3]+2.0)+(x[3]+3.0)^2
# p = 3
# input = [interval(-1,1), interval(0,1), interval(-1,1)]
# Dg = Symbolics.jacobian([g], [x[i] for i=1:p])
# Dg_expr = build_function(Dg, [x[i] for i=1:p])
# my_Dg = eval(Dg_expr[1])
# range_Dg = Base.invokelatest(my_Dg,input)

# println(range_Dg)