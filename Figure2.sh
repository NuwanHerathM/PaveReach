#!/bin/bash

delta_values=(0.0 0.125)
first_value=true
for delta in ${delta_values[@]}
do
    if [ "$first_value" = true ]; then
        first=false
    else
        echo  # Print empty line
    fi
    echo "Classification for delta = $delta"
    docker run -v .:/app pavereach examples/exp_5-1_classification.jl $delta
done