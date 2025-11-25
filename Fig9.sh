#!/bin/bash

filename=ex_5-4_circle_collision

eps_x=0.1
eps_p=0.1
echo "$eps_x"
julia ${filename}.jl $eps_x $eps_p -r --save