#!/bin/bash

echo "Classification of ETCS velocity"
docker run -v .:/app pavereach examples/exp_5-2_classification.jl
