import Dictionaries: gettoken, gettokenvalue, iteratetoken, iteratetoken_reverse, istokenizable, set!

struct FESet{T} <:AbstractIndices{T}
    set::Indices{T}    
end
FESet(v::AbstractVector) = FESet{eltype(v)}(Indices(v))
FESet{T}() where T = FESet(Indices{T}())


function Base.getindex(s::FESet{T},v::SVector{N,I}) where {T,N,I<:Integer} 
    isin,(_,i) = gettoken(s,v)
    if isin
        return gettokenvalue(s,i)
    else 
        isin,(_,i) = gettoken(s,reverse(v))
        if isin
            return gettokenvalue(s,i)
        else
            error("Error: element $v does not belong to the `FESet`")
        end
    end
end

function Base.getindex(s::FESet{T},vector::AbstractVector{SVector{N,I}}) where {T,N,I<:Integer}
    [getindex(s,el) for el in vector]
end

function Base.getindex(s::FESet{T},gen::Base.Generator) where {T}
    (getindex(s,el) for el in gen)
end


# function gettoken(s::FESet{T},v::SVector{N,I}) where {T,N,I<:Integer}             
#     token = gettoken(s.set,v)
#     if !token[1]
#         token = gettoken(s.set,reverse(v))
#     end
#     return token
# end

@inline gettoken(s::FESet{T},v::SVector{N,I}) where {T,N,I<:Integer} = gettoken(s.set,v)
@inline Base.length(s::FESet{T}) where T      = length(s.set)
@inline Base.in(v::SVector{N,I},s::FESet{T}) where {N,T,I<:Integer} = gettoken(s,v)[1] || gettoken(s,reverse(v))[1]
@inline Base.iterate(s::FESet{T}) where T     = iterate(s.set)
@inline Base.iterate(s::FESet{T},i) where T   = iterate(s.set,i)
@inline Base.broadcast(f,s::FESet{T}) where T = broadcast(f,s.set)
@inline iteratetoken(s::FESet{T}) where T     = iteratetoken(s.set)
@inline iteratetoken(s::FESet{T},i) where T   = iteratetoken(s.set,i)
@inline gettokenvalue(s::FESet{T},i) where T        = gettokenvalue(s.set,i)
@inline iteratetoken_reverse(s::FESet{T}) where T   = iteratetoken_reverse(s.set)
@inline iteratetoken_reverse(s::FESet{T},i) where T = iteratetoken_reverse(s.set,i)
@inline istokenizable(s::FESet{T}) where T          = true
@inline add!(s::FESet{T},element::T) where T        = set!(s.set,element)
@inline remove!(s::FESet{T},element::T) where T     = unset!(s.set,element)


