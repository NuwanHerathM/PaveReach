#!/bin/bash

echo "Clear of conflict distance (coarser paving)"
julia exp_5-3_acasxu_clear_distance_coarse.jl

echo # Print empty line
echo "Clear of conflict with velocity change (faster intruder)"
julia exp_5-3_acasxu_safely_clear_distance_faster.jl
