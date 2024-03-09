"""
    HPTuple{N,I<:Integer} 

Abstract supertype for triangles, edges and other sets used in the context of `FEMhp.jl`. 

Every `HPTuple` must implement the function `vals(t::HPTuple)`, returning a vector. An `HPTuple` can be used for indexing: it index according to the output of `vals`. However, for the purposes of hashing and equality checking (vía `isequal` or `==`), it behaves like a set, in the mathematical sense. 

vals can be stored as separated fields or as an array. Additional fields (not returned by `vals`) can be included in an HPTuple type, though it is not recommended. 

# Example: 

    struct MyHPTupleType{I} <: HPTuple{3,I}
        v::SVector{3,I}
        some_property::I
    end
    MyHPTupleType(v::Vector) = MyHPTupleType(SVector{3}(v))
    vals(t::MyHPTupleType) = t.v

    ```julia
    julia> a = MyHPTupleType([1,3,5]);
    julia> b = MyHPTupleType([5,3,1]);
    julia> z = rand(1:10,2,5)
    2×5 Matrix{Int64}:
     6  10  3  10  4
     1   9  5   2  4
    julia> z[:,a]
    2×3 Matrix{Int64}:
     6  3  4
     1  5  4   
    julia> a==b
    true
```
"""
abstract type HPTuple{N,I<:Integer} end

Base.length(t::HPTuple) = length(vals(t))
Base.getindex(t::HPTuple,i) = getindex(vals(t),i)
Base.getindex(t::HPTuple,::Colon) = vals(t)
Base.iterate(t::HPTuple) = iterate(vals(t))
Base.iterate(t::HPTuple, i) = iterate(vals(t),i)
Base.keys(t::HPTuple) = 1:length(t)
Base.Vector(t::HPTuple) = Vector(vals(t))
const hashs_seed = UInt === UInt64 ? 0x793bac59abf9a1da : 0xdea7f1da
function Base.hash(s::HPTuple, h::UInt)
    hv = hashs_seed
    for x in vals(s)
        hv ⊻= hash(x)
    end
    hash(hash(hv, h),hash(typeof(s)))
end

Base.:(==)(t1::T,t2::T) where T<:HPTuple = length(t1)==length(t2) && t1⊆t2
Base.issubset(t1::T,t2::T) where T<:HPTuple = vals(t1) ⊆ vals(t2)


# Using HPTuples for INDEXING:
Base.getindex(A::AbstractArray,t::HPTuple) = getindex(A,vals(t))
Base.getindex(A::AbstractArray,c::Colon,t::HPTuple) = getindex(A,c,vals(t))
Base.getindex(A::AbstractArray,t::HPTuple,c::Colon) = getindex(A,vals(t),c)
Base.view(A::AbstractArray,t::HPTuple) = view(A,vals(t))
Base.view(A::AbstractArray,c::Colon,t::HPTuple) = view(A,c,vals(t))
Base.view(A::AbstractArray,t::HPTuple,c::Colon) = view(A,vals(t),c)


function Base.sizehint!(d::Dictionary,n::Int)
    sizehint!(d.indices.slots,(1+n÷8)*8)
    sizehint!(d.indices.hashes,n)
    sizehint!(d.indices.values,n)
    sizehint!(d.values,n)
end

