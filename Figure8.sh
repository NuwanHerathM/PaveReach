#!/bin/bash

echo "Clear of conflict distance"
docker run -v .:/app pavereach examples/exp_5-3_acasxu_clear_distance.jl

echo # Print empty line
echo "Clear of conflict with velocity change"
docker run -v .:/app pavereach examples/exp_5-3_acasxu_safely_clear_distance.jl
