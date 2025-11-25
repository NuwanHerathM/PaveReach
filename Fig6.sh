#!/bin/bash

filename=ex_5-1_running_example

o_in=1
o_out=1
eps_x=0.001
eps_p=0.001
echo "$o_in, $o_out, $eps"
echo "Subdividing no paramater"
julia ${filename}.jl $o_in $o_out $eps_x --save
echo "Subdividing either existentially or universally quantified parameters"
julia ${filename}.jl $o_in $o_out $eps_x $eps_p -s --save
echo "Combining the subdivision of existentially and universally quantified parameters"
julia ${filename}.jl $o_in $o_out $eps_x $eps_p -r --save