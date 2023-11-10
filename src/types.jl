
import Base: hash, isequal, ==, getindex

using Triangulate
using LinearAlgebra

"""
    TriangleHP{I<:Integer,P<:Integer}

A triangle of an HPMesh, formed by: 

+ `vertices::Tuple{I,I,I}`: indicating the indices of the vertices that form the triangle.
+ `longest::Tuple{I,I}`: the longest edge. 
+ `refine::P`: a number indicating the refinement to be performed.
    - `0`: No refinement
    - `1`: Green refinement (the longest edge is refined)
    - HMMM... think abouot this!
"""

struct TriangleHP{I<:Integer,P<:Integer}
    vertices::Tuple{I,I,I}
    longest::Tuple{I,I}
    refine::Bool
    order::Tuple{P,Bool}
end
TriangleHP(v::V,t,r,o) where V<:AbstractArray = TriangleHP(Tuple(v),t,r,o)

"""
    EdgeHP{I<:Integer,P<:Integer}

An edge of an HPMesh, formed by:

+ `edge::Tuple{I,I}`: indicating the indices of the vertices that form the edge.
+ `p::P`: degree of the polynomial associated with the edge.
+ `r::Bool`: marked for refinement. 
"""
struct EdgeHP{I<:Integer,P<:Integer}
    edge::Tuple{I,I}
    p::P
    m::P
    r::Bool
    EdgeHP{I,P}(v,p,m,r) where {I<:Integer,P<:Integer} = length(v)!=2 ? error("Edge must have two elements, trying to build an edge with $(length(v)) elements instead") : new(v,p,m,r)
end

EdgeHP(v::V,p,m,r) where V<:AbstractArray = EdgeHP(Tuple(v),p,m,r)
"""
struct MeshHP{F<:AbstractFloat,I<:Integer,P<:Integer}
"""
struct MeshHP{F<:AbstractFloat,I<:Integer,P<:Integer}
    points::Matrix{F}
    triangles::Set{TriangleHP{I,P}}
    edges::Set{EdgeHP{I,P}}    
end

function MeshHP(tri::TriangulateIO)
    (;pointlist,trianglelist,edgelist,edgemarkerlist) = tri
    F = eltype(pointlist)
    I = eltype(trianglelist)
    P = UInt8
    build_edge(v::V,m::I) where {V<:AbstractArray, I<:Integer} = EdgeHP{I,P}(v,1,m,false)
    function build_triangle(v::V,p::Matrix{F}) where {V<:AbstractArray,F<:AbstractFloat}
        maxi = argmax(norm(p[:,v[mod1(i+1,3)]]-p[:,v[i]]) for i in 1:3)
        TriangleHP{I,P}(v,(maxi,mod1(maxi,3)),false,(P.(1),true))
    end
    triangles = Set([build_triangle(v,pointlist) for v in eachcol(trianglelist)])
    edges     = Set([build_edge(e,edgemarkerlist[i]) for (i,e) in enumerate(eachcol(edgelist))])
    MeshHP{F,I,P}(pointlist,triangles,edges)
end
"""
Functions that allow creation of a `Set` of edges admitting membership test with tuples. 

    a = EdgeHP((1,2),1,0,false)
    b = EdgeHP((1,2),3,0,true)
    c = EdgeHP((1,3),5,1,false)
    S = Set([a,b,c])
        Set{EdgeHP{Int64,Int64}} with e elements:
            EdgeHP{Int64,Int64}((1,2),3,true)
            EdgeHP{Int64,Int64}((1,3),5,false)
    (1,3)∈S
        true
    (3,1) in S
        true
    
"""
@inline hash(t::Tuple{I,I}) where I<:Integer = hash(t[1])+hash(t[2])
@inline hash(e::EdgeHP) = hash(e.edge)
@inline isequal(e1::EdgeHP,e2::EdgeHP) = isequal(e1.edge,e2.edge) 

@inline isequal(t::Tuple,e2::EdgeHP) = isequal(t,e2.edge) || isequal(reverse(t),e2.edge)
@inline isequal(e2::EdgeHP,t::Tuple) = isequal(t,e2)


#@inline ==(e1::EdgeHP,e2::EdgeHP) = (e1.edge==e2.edge) 
#@inline ==(e::EdgeHP,t::Tuple{I,I}) where I<:Integer = (t==e.edge) 
#@inline ==(t::Tuple{I,I},e::EdgeHP) where I<:Integer = (t==e.edge) 

@inline getindex(e::EdgeHP,k) = e.edge[mod1(k,2)]
@inline getindex(t::TriangleHP,k) = t.vertices[mod1(k,3)]
