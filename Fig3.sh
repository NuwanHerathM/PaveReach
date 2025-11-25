#!/bin/bash

filename=ex_5-1_running_example

for o_in in 1 2
do
    for o_out in 1 2
    do
        eps_x=0.1
        echo "$o_in, $o_out, $eps_x"
        julia ${filename}.jl $o_in $o_out $eps_x --save
    done
done