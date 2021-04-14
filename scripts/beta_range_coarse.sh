#!/bin/bash

for L in 4 10
do
    for b in $(seq 0.8 0.1 3.0)
    do
        julia main.jl $L -n 100000 --beta $b
    done
done