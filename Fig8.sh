#!/bin/bash

filename=ex_5-5_robot_collision

eps_x=0.01
eps_p=0.01
echo "$eps_x"
docker run -v .:/app globalqe examples/${filename}.jl $eps_x $eps_p -r --with_luxor