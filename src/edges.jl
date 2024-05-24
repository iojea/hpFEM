# EDGES
# struct EdgeHP{I} <: HPTuple{2,I}
#     nodes::SVector{2,I}
#     function EdgeHP(v::SVector{2,I}) where I<:Integer
#         if v[1]==v[2]
#             throw("The vertices of an edge must be different. Repeated indices were passed.")
#         else 
#             new{eltype(v)}(v)
#         end
#     end
# end
# EdgeHP(v) = EdgeHP(SVector{2}(v))
# EdgeHP{I}(v) where I<:Integer = EdgeHP(SVector{2,I}(v))
# EdgeHP(a,b) = EdgeHP(SVector{2}([a,b]))
# EdgeHP{I}(a,b) where I<:Integer = EdgeHP{I}(SVector{2,I}([a,b]))
# @inline vals(t::EdgeHP) = t.nodes
# @inline nodes(t::EdgeHP)  = t.nodes


struct EdgeProperties{P<:Integer,Bool} 
    degree::Ref{P}
    marker::P
    refine::Ref{Bool}
end

EdgeProperties(d::P,m::P,r::Bool) where P<:Integer = EdgeProperties(Ref(d),m,Ref(r))
EdgeProperties(d::Ref{P},m::P,r::Bool) where P<:Integer = EdgeProperties(d,m,Ref(r))
EdgeProperties(d::P,m::P,r::Ref{Bool}) where P<:Integer = EdgeProperties(Ref(d),m,r)
EdgeProperties(e::EdgeProperties)  = EdgeProperties(e.degree,e.marker,e.refine)
EdgeProperties{P,Bool}(e::EdgeProperties) where P<:Integer = EdgeProperties(P(e.degree[]),P(e.marker),e.refine)

@inline ismarked(e::EdgeProperties)        = e.refine[]
@inline degree(e::EdgeProperties)          = e.degree[]
@inline marker(e::EdgeProperties)          = e.marker
@inline mark!(e::EdgeProperties)           = e.refine[] = true
@inline setdegree!(e::EdgeProperties,deg) = e.degree[] = deg


            
