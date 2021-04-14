# swendsen wang
# wolff
# single spin flip

using DataStructures
using Statistics
using DelimitedFiles
using Random

include("mc_state.jl")

function ΔE(mc_state::MCState, site::Int, dtheta::Float64)
    # energy difference by rotating a site by dtheta
    NNrotors = mc_state.angle_config[mc_state.neighbours[site]]
    base_angle = mc_state.angle_config[site]
    dE = 0.
    for angle in NNrotors
        dE += -cos(base_angle + dtheta - angle) + cos(base_angle - angle)
    end
    return dE
end

function single_rotation_update!(H::NNXY, mc_state::MCState, dtheta::Float64, β::Float64)
    N = nspins(H)
    angle_config = mc_state.angle_config
    # sweep over all sites
    shuffled_sites = collect(1:N)[randperm(N)]
    for site in shuffled_sites
        dtheta *= rand((-1,1)) # randomly choose to rotate CW or CCW
        Ediff = ΔE(mc_state, site, dtheta)

        if rand() < min(1, exp(-β*Ediff))
            angle_config[site] += dtheta
            angle_config[site] = mod(angle_config[site], 2π)
        end
    end
end

upper_neighbour_periodic(dims::NTuple{2,Int}, site::Int) = site + dims[1] > prod(dims) ? mod(site + dims[1], prod(dims)) : site + dims[1]
#upper_neighbour(dims::NTuple{2,Int}, site::Int) = site + dims[1]

lower_neighbour_periodic(dims::NTuple{2,Int}, site::Int) = site <= dims[1] ? site + prod(dims) - dims[1] : site - dims[1]
#lower_neighbour(dims::NTuple{2,Int}, site::Int) = site - dims[1]

right_neighbour_periodic(dims::NTuple{2,Int}, site::Int) = mod(site, dims[1]) == 0 ? site - dims[1] + 1 : site + 1
#right_neighbour(dims::NTuple{2,Int}, site::Int) = site + 1

left_neighbour_periodic(dims::NTuple{2,Int}, site::Int) = mod(site-1, dims[1]) == 0 ? site + dims[1] - 1 : site - 1
#left_neighbour(dims::NTuple{2,Int}, site::Int) = site - 1