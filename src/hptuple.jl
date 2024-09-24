"""
    TupleHP{L,I}(x::NTuple{L}) 
    TupleHP{L,I}(x1,x2,x3,...)

Construct a statically-sized array `TupleHP`. The type is immutable, so the data must be provided upon construction. The `I` parameter is a type, `I<:Integer`. `L` can be inferred from the data:

    TupleHP{I}(x::Array)

If `I` is not provided the type of the data is inferred from `eltype(x)`:
    TupleHP(x::Array)

`TupleHP` is a subtype of `StaticArray`. As a consequence, an `TupleHP` can be used for indexing, slicing and Linear Algebra operations. However, and `TupleHP` behaves like a mathematical set for the purposes of hashing and equality checking via `isequal`. The goal of this implementation is to use `TupleHP`s as keys of a `Dictionary` allowing search with permutations of the `TupleHP`. 

# Example
```julia
    julia> v = TupleHP(1,2);
    julia> w = TupleHP(2,1);
    julia> isequal(v,w)
    true
    julia> v==w
    false
    julia> using Dictionaries
    julia> d = Dictionary([v],[8])
    1-element Dictionary{TupleHP{Int64},Int64}
    [1, 2] | 8
    julia> w in keys(d)
    true    
    julia> d[w]
    8
```
"""
struct TupleHP{L,I<:Integer} <: StaticArray{Tuple{L},I,1}
    data::NTuple{L,I}

    function TupleHP{L,I}(x::NTuple{L,I}) where {L,I<:Integer}
        check_array_parameters(Tuple{L}, I, Val{1}, Val{L})
        new{L,I}(x)
    end

    function TupleHP{L,I}(x::NTuple{L,Any}) where {L,I}
        check_array_parameters(Tuple{L}, I, Val{1}, Val{L})
        new{L,I}(convert_ntuple(I, x))
    end
end


@inline TupleHP(x::Tuple) = TupleHP{length(x),eltype(x)}(x)
@inline TupleHP{I}(x::Tuple) where {I<:Integer} = TupleHP(convert_ntuple(I, x))

@propagate_inbounds function getindex(v::TupleHP, i::Int)
    getfield(v, :data)[i]
end


const hashs_seed = UInt === UInt64 ? 0x793bac59abf9a1da : 0xdea7f1da
function Base.hash(s::TupleHP, h::UInt)
    hv = hashs_seed
    for x in getfield(s, :data)
        hv ⊻= hash(x)
    end
    hash(hash(hv, h), hash(typeof(s)))
end

Base.isequal(t1::T, t2::T) where {T<:TupleHP} = length(t1) == length(t2) && issubset(t1, t2)


####################
struct DegTuple{I<:Integer} <: StaticArray{Tuple{3},I,1}
    data::NTuple{3,I}

    # function TupleHP{I}(x::NTuple{3,I}) where {I<:Integer}
    #     check_array_parameters(Tuple{3},I,Val{1},Val{3})
    #     new{I}(x)
    # end

    # function TupleHP{I}(x::NTuple{3,Any}) where {I}
    #     check_array_parameters(Tuple{3}, I, Val{1}, Val{3})
    #     new{I}(convert_ntuple(I, x))
    # end
end


# @inline DegTuple{I}(x::Tuple) where I<:Integer = DegTuple(convert_ntuple(I,x))
#@inline DegTuple(x::Tuple) = DegTuple{eltype(x)}(x)

@propagate_inbounds function getindex(v::DegTuple, i::Int)
    getfield(v, :data)[i]
end


function Base.hash(s::DegTuple, h::UInt)
    hv = hashs_seed
    for x in getfield(s, :data)
        hv ⊻= hash(x)
    end
    hash(hash(hv, h), hash(typeof(s)))
end

# Base.:(==)(t1::T,t2::T) where T<:TupleHP = length(t1)==length(t2) && t1⊆t2
Base.isequal(t1::T, t2::T) where {T<:DegTuple} =
    length(t1) == length(t2) && issubset(t1, t2)

# Base.issubset(t1::T,t2::T) where T<:TupleHP = vals(t1) ⊆ vals(t2)

# abstract type HPTuple{N,I<:Integer} end

# Base.length(t::HPTuple) = length(vals(t))
# Base.getindex(t::HPTuple,i) = getindex(vals(t),i)
# Base.getindex(t::HPTuple,::Colon) = vals(t)
# Base.iterate(t::HPTuple) = iterate(vals(t))
# Base.iterate(t::HPTuple, i) = iterate(vals(t),i)
# Base.keys(t::HPTuple) = 1:length(t)
# Base.Vector(t::HPTuple) = Vector(vals(t))
# const hashs_seed = UInt === UInt64 ? 0x793bac59abf9a1da : 0xdea7f1da
# function Base.hash(s::HPTuple, h::UInt)
#     hv = hashs_seed
#     for x in vals(s)
#         hv ⊻= hash(x)
#     end
#     hash(hash(hv, h),hash(typeof(s)))
# end

# Base.:(==)(t1::T,t2::T) where T<:HPTuple = length(t1)==length(t2) && t1⊆t2
# Base.issubset(t1::T,t2::T) where T<:HPTuple = vals(t1) ⊆ vals(t2)


# # Using HPTuples for INDEXING:
# Base.getindex(A::AbstractArray,t::HPTuple) = getindex(A,vals(t))
# Base.getindex(A::AbstractArray,c::Colon,t::HPTuple) = getindex(A,c,vals(t))
# Base.getindex(A::AbstractArray,t::HPTuple,c::Colon) = getindex(A,vals(t),c)
# Base.view(A::AbstractArray,t::HPTuple) = view(A,vals(t))
# Base.view(A::AbstractArray,c::Colon,t::HPTuple) = view(A,c,vals(t))
# Base.view(A::AbstractArray,t::HPTuple,c::Colon) = view(A,vals(t),c)


function Base.sizehint!(d::Dictionary, n::Int)
    sizehint!(d.indices.slots, (1 + n ÷ 8) * 8)
    sizehint!(d.indices.hashes, n)
    sizehint!(d.indices.values, n)
    sizehint!(d.values, n)
end
