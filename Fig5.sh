#!/bin/bash

filename=ex_5-2_stability_controller

eps_x=0.1
echo "$eps_x"
docker run -v .:/app globalqe examples/${filename}.jl $eps_x --with_luxor