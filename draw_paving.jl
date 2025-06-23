using IntervalArithmetic
using Plots; pythonplot()
xlims!((0,4))
ylims!((0,4))

rectangle(p, q) = Shape([p[1],q[1],q[1],p[1]], [p[2],p[2],q[2],q[2]])


plot_1 = plot(rectangle((0,0), (2,2)), color=:blue, legend=:false)
plot_2 = plot(rectangle((2,0), (4,2)), color=:red, legend=:false)

plot(plot_1, plot_2)

function subdivide(cells)
    new_cells = []
    for cell in cells
        boxes = bisect(cell, 0.5)
        push!(new_cells, boxes[1])
        push!(new_cells, boxes[2])
    end
    return new_cells
end

function refine!(grids)
    new_cells = subdivide(grids[end])
    push!(grids, new_cells)
end

cell = IntervalBox(interval(0,4), interval(0,4))
cells = [cell]

grids = [cells]
for i = 1:8
    refine!(grids)
end

target = interval(3,5)

plots = []
for grid in grids
    current_plot = plot()
    for cell in grid
        x = cell[1]
        y = cell[2]
        p = (x.lo, y.lo)
        q = (x.hi, y.hi)
        z = interval(y.lo + x.hi, y.hi+ x.lo)
        if z ⊂ target
            plot!(rectangle(p, q), color=:green1)
        elseif isempty(z ∩ target)
            plot!(rectangle(p, q), color=:cyan)
        else
            plot!(rectangle(p, q), color=:yellow)
        end
    end
    push!(plots, current_plot)
end

plot(plots..., aspect_ratio=1, legend=:false)

gui()