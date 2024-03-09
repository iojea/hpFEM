# TRIANGLES
struct TriangleHP{I} <: HPTuple{3,I}
    nodes::SVector{3,I}
    function TriangleHP(v::SVector{3,I}) where I<:Integer
        if v[1]==v[2] || v[2] == v[3] || v[1]==v[3]
            throw("The vertices of a triangle must be different. Repeated indices were passed.")
        else 
            new{eltype(v)}(v)
        end
    end
end
TriangleHP(v) = TriangleHP(SVector{3}(v))
TriangleHP{I}(v) where I<:Integer = TriangleHP(SVector{3,I}(v))
TriangleHP(a,b,c) = TriangleHP(SVector{3}([a,b,c]))
TriangleHP{I}(a,b,c) where I<:Integer = TriangleHP{I}(SVector{3,I}([a,b,c]))
@inline _eval(t::TriangleHP,k) = t[mod1(k,3)]
@inline edges(t::TriangleHP)   = [EdgeHP(_eval(t,i),_eval(t,i+1)) for i in 1:3]
@inline longestedge(t::TriangleHP) = EdgeHP(_eval(t,1),_eval(t,2))
@inline vals(t::TriangleHP{I}) where I = t.nodes
@inline nodes(t::TriangleHP)  = t.nodes

#const TriangleHP{I} = SVector{3,I} where I<:Integer


# Properties
struct TriangleProperties{P}
    refine::Ref{P}
end
TriangleProperties(val::P) where P<:Integer = TriangleProperties(Ref(val))
TriangleProperties(tp::TriangleProperties) = TriangleProperties(tp.refine)
TriangleProperties{P}(tp::TriangleProperties) where P<:Integer = TriangleProperties(P(tp.refine[]))

@inline ismarked(t::TriangleProperties)    = t.refine[] > 0
@inline isgreen(t::TriangleProperties)     = t.refine[] == 1
@inline isblue(t::TriangleProperties)      = t.refine[] == 2
@inline isred(t::TriangleProperties)       = t.refine[] == 3
@inline mark!(t::TriangleProperties,k)     = t.refine[] = k
@inline mark!(t::TriangleProperties)       = mark!(t,3)


function triangle(t::T,p::AbstractMatrix) where {T<:Union{Tuple,AbstractArray}}
    maxi = argmax(norm(p[:,t[mod1(i+1,3)]]-p[:,t[i]]) for i in 1:3)
    TriangleHP(t[mod1.(maxi:maxi+2,3)])
end

