PATH="../src/"

include(PATH*"hamiltonian.jl")
include(PATH*"mc_state.jl")
include(PATH*"updates.jl")
include(PATH*"measure.jl")

using DrWatson: savename

using ArgParse
using Random
using DelimitedFiles
using JLD2
using JSON
using FileIO

SCRATCH_PATH = "../examples/data/"

global TC_EXACT = 0.89294

function init_mc(parsed_args)
    L = parsed_args["L"]
    beta = parsed_args["beta"]
    seed = parsed_args["seed"]

    t = (beta^(-1) - TC_EXACT)/TC_EXACT

    #=
    if -1 <= t <= -0.5
        dtheta = 0.1
    elseif -0.5 < t <= 0.5
        dtheta = 0.2
    elseif 0.5 < t <= 1.0
        dtheta = 0.5
    else
        dtheta = 0.8
    end 
    
    if -1 <= t <= -0.5
        dtheta = 0.05
    elseif -0.5 < t <= 0.5
        dtheta = 0.1
    elseif 0.5 < t <= 1.0
        dtheta = 0.2
    else
        dtheta = 0.5
    end 
    =#

    Random.seed!(seed)
    H = NNXY((L, L)) # easier to use PBC for calculating spin stiffness
    mc_state = MCState(H)

    MCS = parsed_args["measurements"]
    EQ_MCS = div(MCS, 10)
    skip = parsed_args["skip"]
    
    #mc_opts = (MCS, EQ_MCS, skip, beta, dtheta)
    mc_opts = (MCS, EQ_MCS, skip, beta)

    #d = (L=L, dtheta=dtheta, beta=beta, seed=seed)
    d = (L=L, beta=beta, seed=seed)
    sname = savename(d; digits = 4)

    return H, mc_state, mc_opts, sname
end


function run(parsed_args)
    H, mc_state, mc_opts, sname = init_mc(parsed_args)
    #MCS, EQ_MCS, skip, beta, dtheta = mc_opts
    MCS, EQ_MCS, skip, beta = mc_opts

    energies = zeros(Float64, MCS)
    energies_sqr = zeros(Float64, MCS)

    ρs = zeros(Float64, MCS)

    #println("Running $(typeof(H)), L=$(H.dims[1]), beta=$beta, dtheta=$dtheta")
    println("Running $(typeof(H)), L=$(H.dims[1]), beta=$beta")
    # equil
    for i in 1:EQ_MCS
        single_rotation_update!(H, mc_state, dtheta, beta)
    end

    for i in 1:MCS
        single_rotation_update!(H, mc_state, dtheta, beta)
        energies[i] = energy(H, mc_state)
        energies_sqr[i] = energies[i]^2
        ρs[i] = spin_stiffness(H, mc_state, beta)

        for _ in 1:skip
            single_rotation_update!(H, mc_state, dtheta, beta)
        end
    end

    E = stats_dict(energies)
    E2 = stats_dict(energies_sqr)
    ρ = stats_dict(ρs)

    obs_estimates = Dict{Symbol, Dict}()
    obs_estimates[:energy] = E
    obs_estimates[:sqr_energy] = E2
    obs_estimates[:spin_stiffness] = ρ

    mkpath(SCRATCH_PATH)
    path = joinpath(SCRATCH_PATH, sname)

    observables_file = path * "_observables.json"

    open(observables_file, "w") do io
        JSON.print(io,
            Dict([k => obs_estimates[k] for k in keys(obs_estimates)]),
            2)
    end
end

s = ArgParseSettings()

@add_arg_table! s begin
    "L"
        help = "2D dimension (i.e. L x L)"
        required = true
        arg_type = Int
    "--seed"
        help = "random seed"
        arg_type = Int
        default = 1234
    "--measurements", "-n"
        help = "Number of samples to record"
        arg_type = Int
        default = 100_000
    "--skip", "-s"
        help = "Number of MC steps to perform between each measurement"
        arg_type = Int
        default = 0
    "--beta"
        help = "Inverse temperature"
        arg_type = Float64
        default = 1.
end

parsed_args = parse_args(ARGS, s)

run(parsed_args)