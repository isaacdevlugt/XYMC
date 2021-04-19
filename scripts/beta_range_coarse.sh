#!/bin/bash

for L in 12 16 20 24
do
    for b in 0.7 0.75 0.8 0.85 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2 2.05 2.1 2.15 2.2 2.5 3.0
    do
        julia main.jl $L -n 1000000 --beta $b
    done
done