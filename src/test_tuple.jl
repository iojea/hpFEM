import StaticArraysCore: check_array_parameters,convert_ntuple
import Base: @propagate_inbounds
import Base: getindex
using StaticArrays


struct TupleHP{I<:Integer,L} <: StaticArray{Tuple{L},I, 1}
    data::NTuple{L,I}

    function TupleHP{I,L}(x::NTuple{L,I}) where {I<:Integer,L}
        check_array_parameters(Tuple{L},I,Val{1},Val{L})
        new{I,L}(x)
    end

    function TupleHP{I,L}(x::NTuple{L,Any}) where {I,L}
        check_array_parameters(Tuple{L}, I, Val{1}, Val{L})
        new{I,L}(convert_ntuple(I, x))
    end
end

@inline TupleHP(x::Tuple) = TupleHP{eltype(x),length(x)}(x)
@inline TupleHP{I}(x::Tuple) where I<:Integer = TupleHP{I,length(x)}(convert_ntuple(I,x))

@propagate_inbounds function getindex(v::TupleHP, i::Int)
    getfield(v,:data)[i]
end

const EdgeHP{I} = TupleHP{I,2}
const TriangleHP{I} = TupleHP{I,3} 
