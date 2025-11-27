#!/bin/bash

filename=ex_5-3_dubins

eps_x=0.01
echo "$eps_x"
echo "x"
julia ${filename}.jl $eps_x -x --save
echo "y"
julia ${filename}.jl $eps_x -y --save
echo "theta"
julia ${filename}.jl $eps_x -t --save