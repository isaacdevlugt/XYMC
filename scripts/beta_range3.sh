#!/bin/bash

for L in 12 16 20 24
do
    for b in $(seq 1.215 0.005 1.35)
    do
        julia main.jl $L -n 1000000 --beta $b
    done
done