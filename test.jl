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
        global X_0 = IntervalBox(interval(-10, 10))
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
        global X_0 = IntervalBox(interval(-10, 10))
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
        global X_0 = IntervalBox(interval(0, 6))
        Z = interval(0, 0)
        global intervals = [interval(2, 8), interval(6, 8), Z]
    end
    4 => begin
        # disk example
        n = 1
        p = 3
        @variables x[1:p]
        f_num = [x[1]^2 + x[2]^2 - x[3]]
        f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
        global qe = QuantifierProblem(f_fun, Df_fun, [(Exists, 3)], p, n)
        global X_0 = IntervalBox(interval(-5, 5), interval(-5, 5))
        Z = interval(0, 16)
        global intervals = [Z]
    end
    5 => begin
        # partial example Jirstrand, 1997 (1/2) -> not working yet
        n = 1
        p = 6
        @variables x[1:p]
        f_num = [x[3]*x[1] + x[4] - x[5]*x[1] - x[6]]
        f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
        global qe = QuantifierProblem(f_fun, Df_fun, [(Forall, 3), (Exists, 4), (Exists, 5), (Exists, 6)], p, n)
        global X_0 = IntervalBox([interval(-1, 1), interval(0, 5)])
        Z = interval(0, 0)
        global intervals = [interval(0, 1), interval(-1, 1), interval(0, 1000), Z]
    end
    6 => begin
        # partial example Jirstrand, 1997 (2/2) -> not working yet
        n = 1
        p = 6
        @variables x[1:p]
        f_num = [(x[3]*(x[2] - 1) + 1)^2 - x[5]*(x[2] - 1) - x[6]]
        f_fun, Df_fun = build_function_f_Df(f_num, x, n, p)
        global qe = QuantifierProblem(f_fun, Df_fun, [(Forall, 3), (Exists, 4), (Exists, 5), (Exists, 6)], p, n)
        global X_0 = IntervalBox([interval(-1, 1), interval(0, 5)])
        Z = interval(0, 0)
        global intervals = [interval(0, 1), interval(-1, 1), interval(0, 1000), Z]
    end
    _ => error("Invalid choice")
end

eps = 0.5
# @btime begin
    # global inn, out, delta = pave(qe, intervals, X_0, eps)
    # box = IntervalBox(intervals)
    # is_in = create_is_in(qe, intervals)
    # is_out = create_is_out(qe, intervals)
    # global p_0 = make_membershipcell_root(box, is_in, is_out)
# p_in_0 = make_in_paving(intervals, qe)
# p_out_0 = make_out_paving(intervals, qe)
p_in_0, p_out_0 = make_pz_11(intervals, qe)
inn, out, delta = pave_11(p_in_0, p_out_0, qe, X_0, eps, is_refined=false)
# end


print_inn_out_delta(inn, out, delta)

# Display the results on a plot if requested
if parsed_args["display"]
    display(X_0, inn, out, delta)
end