#!/usr/local/bin/julia

using BenchmarkTools
using IntervalArithmetic

include("genreach2.jl")

@variables x[1:3]
QE5=([x[1]^2+x[2]^2+2*x[1]*x[2]-20*x[1]-20x[2]+100-x[3]],["forall", 1, "forall", 3, "exists", 2],3,1)
print_QE(QE5)


println("Results for ex5 -order 0")
res_QE5_o0=QEapprox_o0(QE5[1], QE5[2], [["forall", 1, "forall", 3, "exists", 2]], QE5[3], QE5[4], [interval(0,6), interval(2,8), interval(6,8)])
println(res_QE5_o0)
