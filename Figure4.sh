#!/bin/bash

echo "Local robustness for p_1"
docker run -v .:/app pavereach examples/exp_5-1_local_robustness_p1.jl

echo # Print empty line
echo "Local robustness for p_2"
docker run -v .:/app pavereach examples/exp_5-1_local_robustness_p2.jl