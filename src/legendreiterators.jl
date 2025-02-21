using StaticArrays

struct LegendreVectorIterator{N,V<:Number}
    x::SVector{N,V}
    function LegendreVectorIterator{N,V}(x) where {N,V}
        all(-1 .<=x .<=1) || throw(DomainException("x must lie in the interval [-1,1]."))
        N==length(x) || throw("N mus coincide with the length of x.")
        new{N,V}(SVector{N,V}(convert.(V,x)))
    end
end
LegendreVectorIterator(x) = LegendreVectorIterator{length(x),eltype(x)}(x)

Base.eltype(::Type{<:LegendreVectorIterator{V}}) where {V} = V
Base.IteratorSize(::Type{<:LegendreVectorIterator}) = Base.IsInfinite()

@inline function Pl_recursion(ℓ,Plm1, Plm2, x)
  # relation is valid from ℓ = 1
  return @. ((2ℓ-1) * x * Plm1 - (ℓ-1) * Plm2)/ℓ
end

function iterate(liv::LegendreVectorIterator)
    Pl   = one.(liv.x)
    Plm1 = zero.(liv.x)
    return Pl,(1,Pl,Plm1)
end 

function iterate(liv::LegendreVectorIterator,state)
    l, Plm1,Plm2 = state
    Pl = Pl_recursion(l, Plm1,Plm2, liv.x)
    return Pl, (l+1, Pl,Plm1)
end

function leg_pol_vec(x,deg)
    liv = LegendreVectorIterator(x)
    first(Iterators.drop(x,deg))
end
