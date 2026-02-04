#!/bin/bash

filename=ex_5-3_dubins

eps_x=0.01
echo "$eps_x"
echo "x"
docker run -v .:/app globalqe examples/${filename}.jl $eps_x -x --with_plots
echo "y"
docker run -v .:/app globalqe examples/${filename}.jl $eps_x -y --with_plots
echo "theta"
docker run -v .:/app globalqe examples/${filename}.jl $eps_x -t --with_plots