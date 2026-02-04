#!/bin/bash

outfile=Tab5.log
rm -f $outfile

filename=ex_5-1_running_example

o_in=1
o_out=1
for eps_x in 0.1 0.01 0.001 0.0001
do
    for eps_p in 0.1 0.01 0.001 0.0001
    do
        echo "$eps_x, $eps_p"
        res=$(docker run -v .:/app globalqe examples/${filename}.jl $o_in $o_out $eps_x $eps_p -r 2>&1)
        echo $res >> $outfile
    done
done