#!/bin/bash

filename=ex_5-5_robot_collision

eps_x=0.1
eps_p=0.1
echo "Subdividing no paramater"
echo "$eps_x"
docker run -v .:/app globalqe examples/${filename}.jl $eps_x --with_luxor
echo "Subdividing either existentially or universally quantified parameters"
echo "$eps_x, $eps_p"
docker run -v .:/app globalqe examples/${filename}.jl $eps_x $eps_p -s --with_luxor
echo "Combining the subdivision of existentially and universally quantified parameters"
echo "$eps_x, $eps_p"
docker run -v .:/app globalqe examples/${filename}.jl $eps_x $eps_p -r --with_luxor