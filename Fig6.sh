#!/bin/bash

filename=ex_5-3_dubins

eps_x=0.01
echo "$eps_x"
echo "x"
julia ${filename}.jl $eps_x -x --with_plots
echo "y"
julia ${filename}.jl $eps_x -y --with_plots
echo "theta"
julia ${filename}.jl $eps_x -t --with_plots