#!/bin/bash

for L in 12 16
do
    for b in $(seq 1.105 0.005 1.21)
    do
        julia main.jl $L -n 1000000 --beta $b
    done
done
