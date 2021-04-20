# Monte Carlo simulations of the 2D classical XY model with periodic boundaries.

To run this code, you must have the following dependencies.

- DrWatson, ArgParse, Random, DelimitedFiles, JSON, FileIO, Measurements, StatsBase

After adding these packages, navigate to `src/` and run `julia main.jl --help` for a list of parameters to specify to run the simulation. Simulations will save a `.json` file with energy, squared energy, and spin stiffness estimates.
