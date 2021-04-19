#!/bin/bash

for L in 12 16 20 24
do
    for b in $(seq 1.005 0.005 1.1)
    do
        julia main.jl $L -n 1000000 --beta $b
    done
done