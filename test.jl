#!/usr/local/bin/julia

using IntervalArithmetic

include("genreach2.jl")
include("pave.jl")
include("quantifierproblem.jl")

using ArgParse

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

# simple example
# n = 1
# p = 3
# @variables x[1:p]
# QE = ([x[1]+x[2]-x[3]],["exists", 2, "exists", 3],p,n)
# Z = interval(3,4)
# intervals = [interval(2,8), Z]

# circle example
n = 1
p = 3
@variables x[1:p]
# QE = ([x[1]^2+x[2]^2-x[3]], ["exists", 2, "exists", 3], p, n)
qe = QuantifierProblem([x[1]^2+x[2]^2-x[3]], [(Exists, 2), (Exists, 3)], p, n)
# println(qe.quantifiers)
# println(negation.(qe.quantifiers))
Z = interval(3, 5)
intervals = [interval(-3, 3), Z]

X_0 = interval(-10, 10)
eps = 0.1
# inn, out, delta = pave(qe, intervals, X_0, eps)
box = IntervalBox(intervals)
is_in = create_is_in(qe, intervals)
is_out = create_is_out(qe, intervals)
p = make_membershipcell_root(box, is_in, is_out)
inn, out, delta = pave(p, qe, X_0, eps)

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