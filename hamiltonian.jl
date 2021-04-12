# 2D Ising PBC

using LinearAlgebra

abstract type AbstractXY end

# TODO: general Ising
# TODO: make sure dims = (x, 1) is interpreted as just dims = x (1D chain)
IntOrTupleInt = Union{Int, NTuple{2,Int}}
BoolOrTupleBool = Union{Int, NTuple{2,Bool}}

struct NNXY <: AbstractXY
    dims::IntOrTupleInt
    PBC::NTuple{2, Bool}
    NNXY(dims,PBC) = new(dims,PBC)
end

nspins(H::AbstractXY) = prod(H.dims)