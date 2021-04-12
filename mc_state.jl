using Distributions

include("hamiltonian.jl")

struct MCState
    angle_config::AbstractVector{Float64}
    neighbours::Array{Array{Int64,1},1}
end

function MCState(H::AbstractXY)
    Ns = nspins(H)
    angle_config = rand(Uniform(0, 2π), Ns)
    neighbours = nearest_neighbours(H.dims, H.PBC)

    return MCState(angle_config, neighbours)
end

function nearest_neighbours(dim::Int, PBC::Bool)
    bonds = [Int[] for _ in 1:dim]
    for i in 2:dim-1
        push!(bonds[i], i+1)
        push!(bonds[i], i-1)
    end

    push!(bonds[1], 2)
    push!(bonds[dim], dim-1)
    if PBC
        push!(bonds[1], dim)
        push!(bonds[dim], 1)
    end

    return bonds
end

function nearest_neighbours(dims::NTuple{2,Int}, PBC::NTuple{2,Bool})
    N = prod(dims)
    bonds = [Int[] for _ in 1:N]
    for i in 1:N
        # left/right, up/low neighbours
        # check if in the outer edges of square
        
        # first row check
        if i <= dims[1]
            # definitely has upper neighbour
            push!(bonds[i], upper_neighbour(dims, i))
            if PBC[2]
                push!(bonds[i], lower_neighbour_periodic(dims, i))
            end

            # first row, last column
            if mod(i, dims[1]) == 0
                # definitely has left neighbour
                push!(bonds[i], left_neighbour(dims, i))
                if PBC[1]
                    push!(bonds[i], right_neighbour_periodic(dims, i))
                end
            
            # first row, first column
            elseif mod(i-1, dims[1]) == 0
                # definitely has right neighbour
                push!(bonds[i], right_neighbour(dims, i))
                if PBC[1]
                    push!(bonds[i], left_neighbour_periodic(dims, i))
                end

            # somewhere in the middle 
            else
                push!(bonds[i], right_neighbour(dims, i))
                push!(bonds[i], left_neighbour(dims, i))
            end

        # last row check
        elseif i > N - dims[1]
            # definitely has lower neighbour
            push!(bonds[i], lower_neighbour(dims, i))

            if PBC[2]
                push!(bonds[i], upper_neighbour_periodic(dims, i))
            end

            # last row, last column
            if mod(i, dims[1]) == 0
                # definitely has left neighbour
                push!(bonds[i], left_neighbour(dims, i))
                if PBC[1]
                    push!(bonds[i], right_neighbour_periodic(dims, i))
                end
            
            # last row, first column
            elseif mod(i-1, dims[1]) == 0
                # definitely has right neighbour
                push!(bonds[i], right_neighbour(dims, i))
                if PBC[1]
                    push!(bonds[i], left_neighbour_periodic(dims, i))
                end

            # somewhere in the middle 
            else
                push!(bonds[i], right_neighbour(dims, i))
                push!(bonds[i], left_neighbour(dims, i))
            end

        # in the bulk
        else
            push!(bonds[i], right_neighbour(dims, i))
            push!(bonds[i], left_neighbour(dims, i))
            push!(bonds[i], lower_neighbour(dims, i))
            push!(bonds[i], upper_neighbour(dims, i))
        end
    end
    return bonds
end