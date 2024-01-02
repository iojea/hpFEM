"""
    EdgeHP{I<:Integer,P<:Integer,Bool}

An edge of an HPMesh, formed by:

+ `nodes::SVector{2,I}`: indicating the indices of the vertices that form the edge. It is implemented as an `SVector` (and not a `Tuple`) in order to be able to use it for indexing, for example, the matrix of points of a given mesh. 
+ `props::EdgeProps{P,Bool}`: Properties associated with the edge.
"""
struct EdgeHP{I<:Integer,P<:Integer,Bool}
    nodes::SVector{2,I}
    degree::Ref{P}
    marker::P
    refine::Ref{Bool}
    #props::EdgeProps{P,Bool}
end

"""
    EdgeHP{I,P}(v::V,p::P,m::P,r::Bool) where {I<:Integer,P<:Integer,V<:AbstractVector}

Construct and `EdgeHP` with vector `v` indicating its `nodes`, degree `p`, marker `m` and refinement indication `r`. `p`, `m` and `r` are packed into a `EdgeProps` struct. 
"""
EdgeHP(v::V,d::P,m::P,r::Bool) where {V<:AbstractArray,P<:Integer} = EdgeHP(SVector{2,eltype(v)}(v),Ref(d),m,Ref(r))
EdgeHP(v::V) where V<:AbstractArray = EdgeHP(v,UInt8(1),UInt8(0),false)
EdgeHP(v::V,d,m,r) where {V<:AbstractArray} = EdgeHP{eltype(v),eltype(promote(d,m)),Bool}(v,d,m,r)
EdgeHP{I,P}(v,d,m,r) where {I<:Integer,P<:Integer} = EdgeHP{I,P,Bool}(v,d,m,r)
@inline ismarked(e::EdgeHP)           = e.refine[]
@inline degree(e::EdgeHP)             = e.degree[]
@inline marker(e::EdgeHP)             = e.marker
@inline nodes(e::EdgeHP)              = e.nodes
@inline mark!(e::EdgeHP)              = e.refine[] = true
@inline set_degree!(e::EdgeHP,deg)    = e.degree[] = deg

#@inline Base.hash(t::SVector{N,I}) where {N,I<:Integer} = sum(Base.hash(node) for node in t)
@inline Base.hash(e::EdgeHP{I,P,Bool}) where {I<:Integer,P<:Integer} = Base.hash(e.nodes)
@inline Base.isequal(e1::EdgeHP,e2::EdgeHP) = Base.isequal(e1.nodes,e2.nodes) 
@inline Base.isequal(t::SVector{2,I},e2::EdgeHP{I,P,Bool}) where {I<:Integer,P<:Integer} = Base.isequal(t,e2.nodes) || Base.isequal(reverse(t),e2.nodes)
@inline Base.isequal(e2::EdgeHP{I,P,Bool},t::SVector{2,I}) where {I<:Integer,P<:Integer} = Base.isequal(t,e2)

# function get_edge(edgeset::Indices{EdgeHP{I,P,Bool}},e::SVector{2,I}) where {I<:Integer,P<:Integer}
#     isin,(_,ind) = gettoken(edgeset,e)
#     if isin
#         return gettokenvalue(edgeset,ind)
#     else
#         isin,(_,ind) = gettoken(edgeset,reverse(e))
#         if isin
#             return gettokenvalue(edgeset,ind)
#         else
#             error("Error: the edge $e does not belong to the edge-set.")
#         end
#     end
# end 


function Base.show(io::IO,e::EdgeHP)
    compact = get(io,:compact,false)
    if !compact
        print(io,"(",e.nodes[1],",",e.nodes[2],")")
        print(io," => (deg = ")
        print(io,degree(e))
        ismarked(e) ? print(io,", refine") : nothing
        if marker(e) == 1
            print(io,", boundary")
        elseif marker(e)>1
            print(io,", interface #$(marker(e)-1)")
        end
        print(io,")")
    else
        show(io,e.nodes)
    end
end

            
