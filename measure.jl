using Measurements
using StatsBase

function stats_dict(f::Function, x::Vector{Float64})
    y = f.(x)
    μ, stdev = mean_and_std(y)
    σ = stdev / sqrt(length(x))
    _, var = mean_and_var(x)
    #μ = mean(y)
    #σ = std(y; mean=μ) / sqrt(length(x))
    dict = Dict("mean" => μ, "error" => σ, "variance" => var)
    return dict
end

function energy(H::NNXY, mc_state::MCState)
    angle_config = mc_state.angle_config
    neighbours = mc_state.neighbours
    Ns = nspins(H)
    z = H.dims isa NTuple ? 4. : 2. # coordination number
    
    E = 0.
    for i in 1:Ns
        θ = angle_config[i]
        θ_NNs = angle_config[neighbours[i]]
        for ϕ in θ_NNs
            E += -cos(θ - ϕ)
        end
    end

    return E / (z*Ns)
end

# TODO: 
function spin_stiffness(H::NNXY, mc_state::MCState, beta::Float64)
    # in the x direction (only taking left/right neighbours)
    # my convention: neighbours[1] = right, neighbours[2] = left
    angle_config = mc_state.angle_config
    neighbours = mc_state.neighbours
    Ns = nspins(H)

    cos_term = 0.
    sin_term = 0.
    for site in 1:Ns
        right, left = neighbours[site][1], neighbours[site][2]
        cos_term += cos(angle_config[site] - right) + cos(angle_config[site] - left)
        sin_term += sin(angle_config[site] - right) + sin(angle_config[site] - left)
    end

    # over-counted pairs in the above loop by 2
    cos_term /= 2
    sin_term /= 2

    return (cos_term - beta*sin_term^2) / Ns
end

stats_dict(x::Vector{Float64}) = stats_dict(identity, x)