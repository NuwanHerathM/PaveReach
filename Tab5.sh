#!/bin/bash

outfile=Tab5.log
rm -f $outfile

filename=ex_5-1_running_example

o_in=1
o_out=1
list_eps_x=( 0.1 0.01 0.001 0.0001 0.00001 )

echo "O^IN_${o_in}, O^OUT_${o_out}"
for eps_x in ${list_eps_x[@]}
do
    echo "$eps_x"
    res=$(julia ${filename}.jl $o_in $o_out $eps_x 2>&1)
    echo $res >> $outfile
done
echo "Results were saved in $outfile."