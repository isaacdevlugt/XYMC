#!/bin/bash

for L in 10 12 16 20 24 30 40 50
do
    for b in $(seq 1.201 0.005 1.35)
    do
        julia main.jl $L -n 1000000 -s 10 --beta $b
    done
done