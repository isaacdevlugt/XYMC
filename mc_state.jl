using Distributions

include("hamiltonian.jl")

struct MCState
    angle_config::AbstractVector{Float64}
    neighbours::Array{Array{Int64,1},1}
end

function MCState(H::AbstractXY)
    Ns = nspins(H)
    angle_config = rand(Uniform(0, 2π), Ns)
    neighbours = nearest_neighbours(H.dims)

    return MCState(angle_config, neighbours)
end

function nearest_neighbours(dim::Int)
    bonds = [Int[] for _ in 1:dim]
    for i in 2:dim-1
        push!(bonds[i], i+1)
        push!(bonds[i], i-1)
    end

    push!(bonds[1], 2)
    push!(bonds[dim], dim-1)
    push!(bonds[1], dim)
    push!(bonds[dim], 1)

    return bonds
end

function nearest_neighbours(dims::NTuple{2,Int})
    # order: upper, lower, right, left
    N = prod(dims)
    bonds = [Int[] for _ in 1:N]
    for i in 1:N
        push!(bonds[i], upper_neighbour_periodic(dims, i))
        push!(bonds[i], lower_neighbour_periodic(dims, i))
        push!(bonds[i], right_neighbour_periodic(dims, i))
        push!(bonds[i], left_neighbour_periodic(dims, i))
    end
    return bonds
end