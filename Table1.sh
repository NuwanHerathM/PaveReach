#!/bin/bash

for d in $(seq 1 4)
do
    echo "Property phi_1 restricted to $d dimension(s)"
    docker run -v .:/app pavereach examples/exp_5-3_acasxu_phi_1_${d}D.jl
done
