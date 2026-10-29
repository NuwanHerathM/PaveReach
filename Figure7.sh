#!/bin/bash

delta_values=(0.98 1.0)

for delta in ${delta_values[@]}
do
    echo "Global robustness for delta = $delta"
    docker run -v .:/app pavereach examples/exp_5-2_global_robustness_landscape.jl $delta
done
