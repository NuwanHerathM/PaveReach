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
    end

    return parse_args(s)
end

parsed_args = parse_commandline()

Base.println("Choose an example to run:")
Base.println("\t1. Simple example")
Base.println("\t\t x_1 | ∃ x_2 ∈ [2, 8], x_1 + x_2 ∈ [3, 4]")
Base.println("\t2. Circle example")
Base.println("\t\t x_1 | ∃ x_2 ∈ [-3, 3], x_1^2 + x_2^2 ∈ [3, 5]")
Base.println("\t3. Jaulin example")
Base.println("\t\t x_1 | ∀ x_3 ∈ [6, 8], ∃ x_2 ∈ [2, 8], x_1^2 + x_2^2 + 2*x_1*x_2 - 20*x_1 - 20*x_2 + 100 - x_3 ∈ [0, 0]")
Base.print("Enter your choice: ")
choice = readline()

# global qe, intervals
@match choice begin
    "1" => begin
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
    "2" => begin
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
    "3" => begin
        # Jaulin example
        n = 1
        p = 4
        @variables x[1:p]
        f_num = [x[1]^2+x[2]^2+2*x[1]*x[2]-20*x[1]-20x[2]+100-x[3]-x[4]]
        f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
        global qe = QuantifierProblem(f_fun, Df_fun, [(Forall, 3), (Exists, 2), (Exists, 4)], p, n)
        global X_0 = interval(0, 6)
        # global X_0 = interval(0.8, 0.9)
        # global X_0 = interval(0.4, 0.45)
        Z = interval(0, 0)
        global intervals = [interval(2, 8), interval(6, 8), Z]
    end
    _ => error("Invalid choice")
end


eps = 0.1
@time begin
    # inn, out, delta = pave(qe, intervals, X_0, eps)
    box = IntervalBox(intervals)
    is_in = create_is_in(qe, intervals)
    is_out = create_is_out(qe, intervals)
    p = make_membershipcell_root(box, is_in, is_out)
    inn, out, delta = pave(p, qe, X_0, eps)
end

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