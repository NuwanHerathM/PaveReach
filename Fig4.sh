#!/bin/bash

filename=ex_5-1_running_example

o_in=1
o_out=1
for eps in 0.1 0.01 0.001
do
    echo "$o_in, $o_out, $eps"
    julia ${filename}.jl $o_in $o_out $eps --save
done