const EdgeHP{I} = TupleHP{2,I}

"""
    EdgeHP{I}(x::Tuple)
    EdgeHP{I}(x1,x2)

constructs an edge formed by a pair of indices of type `I<:Integer`. The type can be inferred from the data:
    
    EdgeHP(x::Tuple)
    EdgeHP(x1,x2)

`EdgeHP{I}` is just an alias for `TupleHP{2,I}`.

"""
EdgeHP(x) = EdgeHP{eltype(x)}(x)

function Base.in(e::T,v::Vector{T}) where T<:EdgeHP
    isin = false
    i    = 1 
    @inbounds while !isin && i≤length(v)
        isequal(e,v[i]) ? isin=true : i+=1
    end
    isin
end

"""
    EdgeProperties(degree::P,marker::P,refine::Bool) where P<:Integer

constructs a `struct` for storing attributes of an edge. These attributes are:
+ `degree`: degree of the polynomial approximator on the edge.
+ `marker`: a marker indicating if the edge belongs to the boundary of the domain, or to an interface or to the interior. 
+ `refine`: `true` if the edge is marked for refinement. 
"""
struct EdgeProperties{P<:Integer,Bool} 
    degree::Ref{P}
    marker::P
    refine::Ref{Bool}
    #adjacent::SVector{2,I}
end

EdgeProperties(d::P,m::P,r::Bool) where P<:Integer = EdgeProperties(Ref(d),m,Ref(r))#,SVector{2,eltype(t)}(t))
EdgeProperties{P}(d,m,r) where P<:Integer = EdgeProperties(Ref(P(d)),P(m),Ref(r))#,SVector{2,eltype(t)}(t))
EdgeProperties(d::Ref{P},m::P,r::Bool) where P<:Integer = EdgeProperties(d,m,Ref(r))#,SVector{2,eltype(t)}(t))
EdgeProperties(d::P,m::P,r::Ref{Bool}) where P<:Integer = EdgeProperties(Ref(d),m,r)#,t)
EdgeProperties(e::EdgeProperties)  = EdgeProperties(e.degree,e.marker,e.refine)#,e.adjacent)
EdgeProperties{P,Bool}(e::EdgeProperties) where P<:Integer = EdgeProperties(P(e.degree[]),P(e.marker),e.refine)#,e.adjacent)


@inline ismarked(e::EdgeProperties)        = e.refine[]
@inline degree(e::EdgeProperties)          = e.degree[]
@inline marker(e::EdgeProperties)          = e.marker
@inline mark!(e::EdgeProperties)           = e.refine[] = true
@inline setdegree!(e::EdgeProperties,deg)  = e.degree[] = deg
@inline isinterior(e::EdgeProperties)      = e.marker != 1
@inline adjacent(e::EdgeProperties)        = e.adjacent


            
