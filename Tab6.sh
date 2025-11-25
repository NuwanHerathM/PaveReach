#!/bin/bash

outfile=Tab6.log
rm -f $outfile

filename=ex_5-1_running_example

o_in=1
o_out=1
list_eps_x=( 0.1 0.01 0.001 0.0001 )

echo "O^IN_${o_in}, O^OUT_${o_out}"
for idx in "${!list_eps_x[@]}"; do
    eps_x=${list_eps_x[$idx]}
    eps_p=$eps_x
    echo "$eps_x, $eps_p"
    echo "Subdividing no paramater"
    res=$(julia ${filename}.jl $o_in $o_out $eps_x $eps_p 2>&1)
    echo $res >> $outfile
    echo "Subdividing either existentially or universally quantified parameters"
    res=$(julia ${filename}.jl $o_in $o_out $eps_x $eps_p -s 2>&1)
    echo $res >> $outfile
    echo "Combining the subdivision of existentially and universally quantified parameters"
    res=$(julia ${filename}.jl $o_in $o_out $eps_x $eps_p -r 2>&1)
    echo $res >> $outfile
done
echo "Results were saved in $outfile."