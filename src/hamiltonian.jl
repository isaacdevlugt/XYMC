# 2D XY Model PBC
# TODO: OBC

abstract type AbstractXY end

IntOrTupleInt = Union{Int, NTuple{2,Int}}
#BoolOrTupleBool = Union{Int, NTuple{2,Bool}}

struct NNXY <: AbstractXY
    dims::IntOrTupleInt
    #PBC::NTuple{2, Bool}
    #NNXY(dims,PBC) = new(dims,PBC)
    NNXY(dims) = new(dims)
end

nspins(H::AbstractXY) = prod(H.dims)