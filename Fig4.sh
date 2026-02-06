#!/bin/bash

filename=ex_5-1_running_example

o_in=1
o_out=1
eps_x=0.001
eps_p=0.001
echo "$o_in, $o_out, $eps"
echo "Subdividing no paramater"
docker run -v .:/app globalqe examples/${filename}.jl $o_in $o_out $eps_x --with_luxor
echo "Subdividing either existentially or universally quantified parameters"
docker run -v .:/app globalqe examples/${filename}.jl $o_in $o_out $eps_x $eps_p -s --with_luxor
echo "Combining the subdivision of existentially and universally quantified parameters"
docker run -v .:/app globalqe examples/${filename}.jl $o_in $o_out $eps_x $eps_p -r --with_luxor