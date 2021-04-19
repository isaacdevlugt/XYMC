#!/bin/bash

for L in 20 24
do
    for b in $(seq 0.9 0.005 1)
    do
        julia main.jl $L -n 1000000 --beta $b
    done
done
