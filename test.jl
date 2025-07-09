#!/usr/local/bin/julia

using IntervalArithmetic

include("genreach2.jl")
include("pave.jl")
include("quantifierproblem.jl")

using ArgParse

using Match

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--display", "-d"
            help = "display the results in a plot"
            action = :store_true
        "--example"
            help = "choose an example to run"
            arg_type = Int
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

if !isnothing(parsed_args["example"])
    global choice = parsed_args["example"]
else
    Base.println("Choose an example to run:")
    Base.println("\t1. Simple example")
    Base.println("\t\t x_1 | ∃ x_2 ∈ [2, 8], x_1 + x_2 ∈ [3, 4]")
    Base.println("\t2. Circle example")
    Base.println("\t\t x_1 | ∃ x_2 ∈ [-3, 3], x_1^2 + x_2^2 ∈ [3, 5]")
    Base.println("\t3. HSVJ05 example")
    Base.println("\t\t x_1 | ∀ x_3 ∈ [6, 8], ∃ x_2 ∈ [2, 8], x_1^2 + x_2^2 + 2*x_1*x_2 - 20*x_1 - 20*x_2 + 100 - x_3 ∈ [0, 0]")
    Base.print("Enter your choice: ")
    global choice = parse(Int, readline())
end

# global qe, intervals
@match choice begin
    1 => begin
        # simple example
        n = 1
        p = 3
        @variables x[1:p]
        f_num = [x[1]+x[2]-x[3]]
        f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
        global qe = QuantifierProblem(f_fun, Df_fun, [(Exists, 2), (Exists, 3)], p, n)
        global X_0 = interval(-10, 10)
        Z = interval(3, 4)
        global intervals = [interval(2, 8), Z]
    end
    2 => begin
        # circle example
        n = 1
        p = 3
        @variables x[1:p]
        f_num = [x[1]^2+x[2]^2-x[3]]
        f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
        global qe = QuantifierProblem(f_fun, Df_fun, [(Exists, 2), (Exists, 3)], p, n)
        global X_0 = interval(-10, 10)
        Z = interval(3, 5)
        global intervals = [interval(-3, 3), Z]
    end
    3 => begin
        # HSVJ05 example (Quantified set inversion algorithm with applications to control, 2005)
        n = 1
        p = 4
        @variables x[1:p]
        f_num = [x[1]^2+x[2]^2+2*x[1]*x[2]-20*x[1]-20x[2]+100-x[3]-x[4]]
        f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
        global qe = QuantifierProblem(f_fun, Df_fun, [(Forall, 3), (Exists, 2), (Exists, 4)], p, n)
        global X_0 = interval(0, 6)
        # global X_0 = interval(3.75, 3.9375)
        # global X_0 = interval(3.75, 4) # is inside
        # global X_0 = interval(3.5, 4) # is not inside but should :(
        Z = interval(0, 0)
        global intervals = [interval(2, 8), interval(6, 8), Z]
    end
    _ => error("Invalid choice")
end

eps = 0.1
@btime begin
    # global inn, out, delta = pave(qe, intervals, X_0, eps)
    # box = IntervalBox(intervals)
    # is_in = create_is_in(qe, intervals)
    # is_out = create_is_out(qe, intervals)
    # global p_0 = make_membershipcell_root(box, is_in, is_out)
    p_0 = make_paving(intervals, qe)
    global inn, out, delta = pave(p_0, qe, X_0, eps)
end




# cell = make_paving(intervals, qe)
# println(is_in(cell, X_0, qe.quantifiers, 1))
# println(is_out(cell, X_0, qe.quantifiers, 1))



# is_in = create_is_in(qe, intervals)
# is_out = create_is_out(qe, intervals)

# pos = [i for (q, i) in qe.quantifiers]
# real_pos = pos .- 1
# invpermute!(intervals, real_pos)

# cell = make_paving(intervals, is_in, is_out)


# intervals = [interval(2, 8),interval(8,9),interval(0, 0)]
# is_in = create_is_in(qe, intervals)
# is_out = create_is_out(qe, intervals)

# pos = [i for (q, i) in qe.quantifiers]
# real_pos = pos .- 1
# invpermute!(intervals, real_pos)
# push!(cell, intervals, is_in, is_out)

# intervals = [interval(2, 8),interval(6,8),interval(0, 0)]
# is_in = create_is_in(qe, intervals)
# is_out = create_is_out(qe, intervals)

# pos = [i for (q, i) in qe.quantifiers]
# real_pos = pos .- 1
# invpermute!(intervals, real_pos)
# push!(cell, intervals, is_in, is_out)
# push!(cell, [interval(6, 8),interval(2, 8),interval(0, 0)])
# push!(cell, [interval(6, 8),interval(2, 3),interval(0, 0)])
# remove!(cell, [interval(6, 8),interval(2, 8),interval(0, 0)])                 # does something
# remove!(cell, [interval(6, 8),interval(2, 8),interval(0, 0), interval(8,9)])  # does nothing as expected
# remove!(cell, [interval(6, 8),interval(2, 8)])                                # does nothing as expected
# print_tree(cell)





# remove!(cell, [interval(8,9),interval(2,5),interval(0, 0)])
# print_tree(cell)
# bisect!(cell, qe)
# bisect!(cell, qe)
# remove!(cell, [interval(8,9),interval(2,5),interval(0, 0)])
# print_tree(cell)



# using Plots; pythonplot()
# cell_0 = inn
# xs = cell_0.box[1]
# ys = cell_0.box[2]
# xlims!((xs.lo, xs.hi))
# ylims!((ys.lo, ys.hi))

# rectangle(p, q) = Shape([p[1],q[1],q[1],p[1]], [p[2],p[2],q[2],q[2]])

# current_plot = plot()

# function draw_mcell(cell)
#     if width(cell) < 4
#         xs = cell.box[1]
#         ys = cell.box[2]
#         p = (xs.lo, ys.lo)
#         q = (xs.hi, ys.hi)

#         plot!(rectangle(p, q), color=cell.is_out(interval(3.75,3.875)) ? :cyan : :yellow)
#     else
#         for child in cell.children
#             draw_mcell(child)
#         end
#     end
# end

# draw_mcell(cell_0)

# plot(current_plot, aspect_ratio=1, legend=:false)
# gui()


# println(height(p_0))
# println(p_0.is_in(X_0))
# println(p_0.is_out(X_0))
# # k = 0
# for k = 1:17
#     # global k += 1    
#     bisect!(p_0, qe)
# end
# println(height(p_0))
# println(width(p_0))
# println(p_0.is_in(X_0))
# println(p_0.is_out(X_0))

# print_tree(p_0)

print_inn_out_delta(inn, out, delta)

# Display the results on a plot if requested
if parsed_args["display"]

    using Plots; pythonplot()
    ylims!((-0.1,0.1))
    xticks!((X_0.lo:1:X_0.hi))

    function draw_lines(intervals, color)
        ys = [0, 0]
        for interval in intervals
            xs = [interval.lo, interval.hi]
            plot!(xs, ys, linewidth=2, color=color, legend=:false)
        end
    end

    draw_delta_lines(delta) = draw_lines(merge_intervals(delta), :yellow)
    draw_inn_lines(inn) = draw_lines(merge_intervals(inn), :green)
    draw_out_lines(out) = draw_lines(merge_intervals(out), :cyan)

    draw_delta_lines(delta)
    draw_inn_lines(inn)
    draw_out_lines(out)
    gui()

end