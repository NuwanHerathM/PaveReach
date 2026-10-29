#!/bin/bash

echo "Global robustness of class 1"
docker run -v .:/app pavereach examples/exp_5-1_global_robustness_landscape_class_1.jl

echo # Print empty line
echo "Global robustness of class 2"
docker run -v .:/app pavereach examples/exp_5-1_global_robustness_landscape_class_2.jl

echo # Print empty line
echo "Global robustness of both classes"
docker run -v .:/app pavereach examples/exp_5-1_global_robustness_landscape_complete.jl
