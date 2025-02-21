using StaticArrays

struct LIVec{V<:Number,N}
    x::SVector{N,V}
    function LIVec{V,N}(x) where {V,N}
        all(-1 .<=x .<=1) || throw(DomainException("x must lie in the interval [-1,1]."))
        N==length(x) || throw("N mus coincide with the length of x.")
        new{V,N}(SVector{N,V}(convert.(V,x)))
    end
end
LIVec(x) = LIVec{eltype(x),length(x)}(x)

@inline function Pl_recursion(ℓ,Plm1, Plm2, x)
  # relation is valid from ℓ = 1
  return @. ((2ℓ-1) * x * Plm1 - (ℓ-1) * Plm2)/ℓ
end

function Base.iterate(liv::LIVec)
    Pl   = one.(liv.x)
    Plm1 = zero.(liv.x)
    return Pl,(1,Pl,Plm1)
end 

function Base.iterate(liv::LIVec,state)
    l, Plm1,Plm2 = state
    Pl = Pl_recursion(l, Plm1,Plm2, liv.x)
    return Pl, (l+1, Pl,Plm1)
end

function leg_pol_vec(x,deg)
    liv = LIVec(x)
    first(Iterators.drop(liv,deg))
end


struct LI{V<:Number}
    x::V
    function LI{V}(x) where {V}
        all(-1 <=x <=1) || throw(DomainException("x must lie in the interval [-1,1]."))
        new{V}(convert(V,x))
    end
end
LI(x) = LI{typeof(x)}(x)


function Base.iterate(li::LI)
    Pl   = one(li.x)
    Plm1 = zero(li.x)
    return Pl,(1,Pl,Plm1)
end 

function Base.iterate(li::LI,state)
    l, Plm1,Plm2 = state
    Pl = Pl_recursion(l, Plm1,Plm2, li.x)
    return Pl, (l+1, Pl,Plm1)
end

function leg_pol(x,deg)
    li = LI(x)
    first(Iterators.drop(li,deg))
end


struct LegendreIterator{V<:Number,N}
    x::MVector{N,V}
    y::MVector{N,V}
    z::MVector{N,V}
    buffer::MVector{N,V}
    function LegendreIterator{V,N}(x) where {V,N}
        all(-1 .<=x.<=1) || throw("x must lie in the interval [-1,1].")
        N==length(x) || throw("N must coincide with the length of x.")
        # n>=0 || throw("n must be an integer n≥1. n=$n was passed.")
        y = MVector{length(x),V}(convert.(V,x))
        new{V,length(x)}(y,similar(y),similar(y),similar(y)) 
    end
end
LegendreIterator(x) = LegendreIterator{eltype(x),length(x)}(x)

Base.IteratorSize(::Type{<:LegendreIterator}) = Base.IsInfinite()


function Base.iterate(iter::LegendreIterator)
    iter.y .= zero(eltype(iter.x))
    iter.z .= one(eltype(iter.x))
    return iter.z,0
end

function Base.iterate(iter::LegendreIterator,k)
    iter.buffer .= iter.y
    iter.y .= iter.z
    iter.z .= Pl_recursion(k+1,iter.z,iter.buffer,iter.x)
    return iter.z,k+1
end


function pol_leg(x,n)
    li = LegendreIterator(x)
    first(Iterators.drop(li,n))
end

struct StandardBasis{T<:Number,I<:Integer}
    x::SVector{2,T}
    p₁::I
    p₂::I
    p₃::I
    function StandardBasis{T,I}(x,p₁,p₂,p₃) where {T,I}
        -1 .<= x .<= 1 || throw(DomainError("x must lie between -1 and 1."))
        new{T,I}(SVector{2,T}(convert.(T,x)),p₁,p₂,p₃)
    end
end

@inline function iterate(basis::StandardBasis)
    p = eltype(x)
    state = (0x00,0x00)
    return p,(p,state)
end

  return @. ((2ℓ-1) * x * Plm1 - (ℓ-1) * Plm2)/ℓ

function iterate(basis::StandardBasis,state)
    i,j = state[2]
    if j<=min(p₂,p₃-i)
        
end
    for i in 0:p₁
        for j in 0:min(p₂,p₃-i)

@inline function iterate(li::Legendre2DIterator, state)
    valid, I = __inc(state.I, iter.indices)
    valid || return nothing
    return CartesianIndex(I...), CartesianIndex(I...)
end


@inline function __inc(state::Tuple{Int}, indices::Tuple{OrdinalRangeInt})
    rng = indices[1]
    I = state[1] + step(rng)
    valid = state[1] != last(rng)
    return valid, (I,)
end
@inline function __inc(state::Tuple{Int,Int,Vararg{Int}}, indices::Tuple{OrdinalRangeInt,OrdinalRangeInt,Vararg{OrdinalRangeInt}})
    rng = indices[1]
    I = state[1] + step(rng)
    if state[1] != last(rng)
        return true, (I, tail(state)...)
    end
    valid, I = __inc(tail(state), tail(indices))
    return valid, (first(rng), I...)
end



