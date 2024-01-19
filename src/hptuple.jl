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
# import Dictionaries: gettokenvalue, settokenvalue!,gettoken,gettoken!, deletetoken!, iteratetoken,distinct,pairs

# struct DictHP{N,I<:Integer,V}
#     dict::AbstractDictionary{SVector{N,I},V}
# end


# DictHP(ind::Indices{SVector{N,I}},vals::V) where {N,I<:Integer,V} = DictHP(Dictionary(ind,vals))

# function DictHP(iterable,vals) 
#     if !(eltype(iterable)<:SVector{N,I} where {N,I<:Integer})
#         el = first(Iterators.take(iterable,1))
#         N  = length(el)
#         I  = eltype(el)
#         try 
#             iterable = SVector{N,I}.(iterable)
#         catch 
#             throw("Elements in iterable cannot be converted to SVector{N,I}.")
#         end
#     end
#     N = length(eltype(iterable))
#     if N == 2 && length(union(iterable,reverse.(iterable))) != 2length(iterable)
#         throw("Repeated elements in list of keys. Take into account that `[a,b]` and `[b,a]` are considered as the same edge.")
#     elseif N==3 && length(union(iterable,circshift.(iterable,1),circshift.(iterable,2))) != 3length(iterable)
#         throw("Repeated elements in list of keys. Take into account that `[a,b,c]`, `[c,a,b]` and `[b,c,a]` are considered as the same triangle.")
#     end
#     DictHP(Indices(iterable),vals)
# end

# DictHP{N,I}(iterable,vals) where {N,I<:Integer} = DictHP{N,I,eltype(vals)}(SVector{N,I}.(iterable),vals)

# DictHP{N,I,V}(iterable,vals) where {N,I<:Integer,V} = DictHP(Indices(SVector{N,I}.(iterable)),V.(vals))

# function DictHP{N,I,V}() where {N,I<:Integer,V} 
#     inds = Vector{SVector{N,I}}()
#     vals = Vector{V}()
#     DictHP(inds,vals)
# end

# function dicthp(iter)
#     DictHP(dictionary(iter))
# end

# @inline Base.keys(dict::DictHP) = keys(dict.dict)
# @inline Base.vals(dict::DictHP) = vals(dict.dict)
# @inline gettokenvalue(dict::DictHP,token) = gettokenvalue(dict.dict,token)
# @inline gettoken(dict::DictHP{N,I,V},el::SVector{N,I}) where {N,I<:Integer,V} = gettoken(dict.dict,el)
# @inline gettoken(dict::DictHP{N,I,V},el) where {N,I<:Integer,V} = gettoken(dict,SVector{N,I}(el))
# @inline gettoken!(dict::DictHP{N,I,V},el::SVector{N,I}) where {N,I<:Integer,V} = gettoken!(dict.dict,el)
# # @inline set!(dict::DictHP{N,I,V},el::SVector{N,I},val::V) where {N,I<:Integer,V} = set!(dict.dict,el,val)
# # @inline unset!(dict::DictHP{N,I,V},el::SVector{N,I}) where {N,I<:Integer,V} = set!(dict.dict,el)
# @inline Base.iterate(dict::DictHP)   = iterate(dict.dict)
# @inline Base.iterate(dict::DictHP,i) = iterate(dict.dict,i) 
# @inline iteratetoken(dict::DictHP)   = iteratetoken(dict.dict)
# @inline iteratetoken(dict::DictHP,i) = iteratetoken(dict.dict,i)
# @inline Base.broadcast(f,dict::DictHP) = DictHP(broadcast(dict.dict,i))
# @inline add!(d::DictHP{N,I,V},el,val) where {N,I<:Integer,V} = set!(d.dict,el,V(val))
# @inline remove!(d::DictHP{N,I,V},el) where {N,I<:Integer,V} = unset!(d.dict,el)

# @inline Base.in(v,dict::DictHP{2,I,V}) where {I<:Integer,V} = gettoken(dict,v)[1] || gettoken(dict,reverse(v))[1]

# @inline Base.in(v,dict::DictHP{3,I,V}) where {I<:Integer,V} = gettoken(dict,v)[1] || gettoken(dict,circshift(v,1))[1] || gettoken(dict,circshift(v,2))[2] 

# @inline Base.map(f,d::DictHP) = DictHP(map(f,d.dict))
# @inline Base.filter(f,d::DictHP) = DictHP(filter(f,d.dict))
# @inline Base.filter!(f,d::DictHP) = filter!(f,d.dict)
# @inline Base.findall(f,d::DictHP) = DictHP(f,d.dict)
# @inline Base.fill(val,d::DictHP) = DictHP(fill(val,d.dict))
# @inline pairs(d::DictHP) = DictHP(pairs(d.dict))

# function Base.getindex(dict::DictHP{2,I,V},v) where {I<:Integer,V}
#     isin,token = gettoken(dict,v)
#     if !isin 
#         isin, token = gettoken(dict,reverse(v))
#     end
#     if isin 
#         return gettokenvalue(dict,token)
#     else
#         throw("edge does not belong to dict")
#     end
# end


# function Base.getindex(dict::DictHP{3,I,V},v) where {I<:Integer,V}
#     isin,token = gettoken(dict,v)
#     if !isin
#         isin,token = gettoken(dict,circshift(v,1))
#         if !isin
#             isin,token = gettoken(dict,circshift(v,2))
#         end
#     end
#     if isin 
#         return gettokenvalue(dict,token)
#     else
#         throw("triangle does not belong to dict")
#     end
# end


# # function Base.sizehint!(ind::Indices,n)
# #     sizehint!(ind.slots,n)
# #     sizehint!(ind.hashes,n)
# #     sizehint!(ind.vals,n)
# # endc# function Base.sizehint!(dict::DictHP,n)
# #     sizehint!(keys(dict),n)
# #     sizehint!(_vals(dict),n)
# # end

# function Base.show(io::IO,dict::DictHP)
#     show(io,dict.dict)
# end

# function Base.show(io::IO, mime::MIME"text/plain", d::DictHP)
#     show(io, mime,d.dict)
# end
