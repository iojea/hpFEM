using Base: _bool

import Base: hash, isequal, ==, getindex, iterate

using Triangulate
using LinearAlgebra
using Printf
using StaticArrays
using Plots
import Plots.plot

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
struct TriangleHP{I<:Integer,P<:Integer,B<:Bool}
    vertices::NTuple{3,I}
    longest::NTuple{2,I}
    refine::B
    order::Tuple{P,B}
    function TriangleHP(v,l,r,o)
        if length(v)!=3
            error("A triangle must have 3 edges, $(length(v)) were passed.")
        elseif !(l ⊆ v)
            error("Longest edge must be contained in the triangle. $l was passed as an edge of triangle $v")
        elseif o[1]∉1:3
            error("`order` must be of type `Tuple{P,Bool} where P<:Integer`. The tuple gives an ordering for reading the triangle such that the degree of the polynomials on the edges is increasing.The number indicates the index of a vertex and the boolean value indicates whether the triangle must be read according to its possitive orientation.")
        else 
            I = eltype(v); P = typeof(o[1])
            new{I,P,Bool}(v,l,r,o)
        end
    end
end


"""
    TriangleHP(v,l,r,o)

Builds a triangle with vertices defined in the `Tuple` or `Array` `v`, longest edge defined in the `Tuple` `l`, refinement indication in the 

"""
TriangleHP(v::V,l,r,o) where V<:AbstractArray = TriangleHP(Tuple(v),l,r,o)
TriangleHP(v::V,l,r::Bool) where V<:Union{AbstractArray,Tuple} = TriangleHP(Tuple(v),l,r,(1,true))
TriangleHP(v::V,l) where V<:Union{AbstractArray,Tuple} = TriangleHP(Tuple(v),l,false)
TriangleHP(v::V,l,o::Tuple{P,Bool}) where {V<:Union{AbstractArray,Tuple},P<:Integer} = TriangleHP(Tuple(v),l,false,o)


"""
    EdgeHP{I<:Integer,P<:Integer}

An edge of an HPMesh, formed by:

+ `edge::Tuple{I,I}`: indicating the indices of the vertices that form the edge.
+ `p::P`: degree of the polynomial associated with the edge.
+ `r::Bool`: marked for refinement. 
"""
struct EdgeHP{I<:Integer,P<:Integer,B<:Bool}
    edge::NTuple{2,I}
    p::P
    m::P
    r::B
    EdgeHP{I,P,B}(v,p,m,r) where {I<:Integer,P<:Integer,B<:Bool} = length(v)!=2 ? error("Edge must have two elements, trying to build an edge with $(length(v)) elements instead") : new(v,p,m,r)
end

"""
    EdgeList{EdgeHP{I<:Integer,P<:Integer,B<:Bool}}

a set for storing edges. Some special methods are implemented for `EdgeList`s in order to simplify access to the edges stored in the set. 
"""
struct EdgeList{I<:Integer,P<:Integer,B<:Bool}
    set::Set{EdgeHP{I,P,B}}
end
"""
struct MeshHP{F<:AbstractFloat,I<:Integer,P<:Integer}
"""
struct MeshHP{F<:AbstractFloat,I<:Integer,P<:Integer,B<:Bool}
    points::Matrix{F}
    triangles::Set{TriangleHP{I,P,B}}
    edges::EdgeList{I,P,B}    
end

function MeshHP(tri::TriangulateIO)
    (;pointlist,trianglelist,edgelist,edgemarkerlist) = tri
    F = eltype(pointlist)
    I = eltype(trianglelist)
    P = UInt8
    build_edge(v::V,m::I) where {V<:AbstractArray, I<:Integer} = EdgeHP{I,P,Bool}(Tuple(v),1,m,false)
    function build_triangle(v::V,p::Matrix{F}) where {V<:AbstractArray,F<:AbstractFloat}
        maxi = argmax(norm(p[:,v[mod1(i+1,3)]]-p[:,v[i]]) for i in 1:3)
        TriangleHP(v,(v[maxi],v[mod1(maxi+1,3)]),false,(P(1),true))
    end
    triangles = Set([build_triangle(v,pointlist) for v in eachcol(trianglelist)])
    edges     = EdgeList(Set([build_edge(e,edgemarkerlist[i]) for (i,e) in enumerate(eachcol(edgelist))]))
    MeshHP(pointlist,triangles,edges)
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

"""
    FEMhp.hash(t::NTuple{N,I}) where {N,I<:Integer} 

Computes the hash of a `Tuple` `t` as the sum of the hashes of the members of `t`. This method is introduced in order to be able to determine if an edge or triangle belongs to a triangulation by using only the tuple as key.
"""
@inline hash(t::NTuple{N,I}) where {N,I<:Integer} = sum(hash(e) for e in t)
@inline hash(e::EdgeHP) = hash(e.edge)
@inline hash(t::TriangleHP) = hash(t.vertices)


@inline isequal(e1::EdgeHP,e2::EdgeHP) = isequal(e1.edge,e2.edge) 
@inline isequal(t::NTuple{2,I},e2::EdgeHP) where I<:Integer = isequal(t,e2.edge) || isequal(reverse(t),e2.edge)
@inline isequal(e2::EdgeHP,t::NTuple{2,I}) where I<:Integer = isequal(t,e2)


@inline isequal(t1::TriangleHP,t2::TriangleHP) = isequal(t1.vertices,t2.vertices) 
@inline isequal(t::NTuple{3,I},tri::TriangleHP) where I<:Integer = isequal(t,tri.vertices) || isequal(reverse(t),tri.vertices)
@inline isequal(tri::TriangleHP,t::NTuple{3,I}) where I<:Integer = isequal(t,tri)

@inline getindex(e::EdgeHP,k) = e.edge[mod1(k,2)]
@inline getindex(t::TriangleHP,k) = t.vertices[mod1(k,3)]

@inline is_marked_for_refinement(t::TriangleHP) = t.m
@inline is_marked_for_refinement(e::EdgeHP) = e.m
@inline degree(e::EdgeHP) = e.p
@inline get_mark(e::EdgeHP) = e.m


iterate(el::EdgeList) = iterate(el.set)
iterate(el::EdgeList,i) = iterate(el.set,i)

read_in_order(t::TriangleHP) = t.order[2] ? (t[i] for i in t.order[1]:t.order[1]+2) : (t[i] for i in t.order[1]:-1:t.order[1]-2)

function get_edges(t::TriangleHP)  
    if t.order[2] 
        return ((t[i],t[i+1]) for i in t.order[1]:t.order[1]+2) 
    else 
        ((t[i],t[i-1]) for i in t.order[1]:-1:t.order[1]-2)
    end
end

function get_edge(t::NTuple{2,I},el::EdgeList{I,P,B}) where {I<:Integer,P<:Integer,B<:Bool}
    getkey(el.set.dict[t])
end

function pop_edge(t::NTuple{2,I},el::EdgeList{I,P,B}) where {I<:Integer,P<:Integer,B<:Bool}
    a = get_edge(t,el)
    delete!(el,t)
    return a
end

function get_degrees(t::TriangleHP{I,P,B},edges::EdgeList{I,P,B}) where {I<:Integer,P<:Integer,B<:Bool}
    (get_edge(t,el).p for t in get_edges(t))
end


function circular_mesh(c,r,h)
    n      = Int(1+2π*r÷h)
    θ      = range(start=0,stop=2π,length=n)[1:end-1]
    points = hcat([[c[1]+cos(t),c[2]+sin(t)] for t in θ]...)
    edges  = hcat([[i,i+1] for i in 1:length(points)-1]...,
                          [length(points),1])
    marks  = ones(length(edges))
    tri = TriangulateIO(;pointlist=points,edgelist=edges,edgemarkerlist=marks)
    maxa = Printf.@sprintf "%.15f" h^2/2
    angle= Printf.@sprintf "%.15f" 30.
    (tri,) = triangulate("ea$(maxa)q$(angle)Q",tri)
    mesh = MeshHP(tri)
end;
circular_mesh(r,h) = circular_mesh(zeros(2),r,h)
circular_mesh(h)   = circular_mesh(1,h)


"""

    plot(mesh::HPMesh,color_scheme::Symbol)

plots `mesh`. Border edges are thicker. Color indicates the degree asigned to the edge.
"""
function plot(mesh::MeshHP;color_scheme::Symbol=:blues,cb=true,kwargs...)
    (;points,edges) = mesh
    nc      = Int(maximum(e.p for e in edges))
    mc      = Int(minimum(e.p for e in edges))
    pal     = palette(color_scheme,nc-mc+2)#nc-mc>0 ? palette(color_scheme,nc-mc+1) : palette(color_scheme,1)
    plt     = plot(bg=:gray30,legend=:none,border=:none;kwargs...)
    l = @layout [a{0.95w} b]
    F = eltype(mesh.points)
    x = MMatrix{2,2}(zeros(F,2,2))
    for e in edges.set
        x .= [points[:,e[1]] points[:,e[2]]]
        if e.m==1
            plot!(plt,x[1,:],x[2,:],c=pal[e.p-mc+1],linewidth=4)
        else
            plot!(plt,x[1,:],x[2,:],c=pal[e.p-mc+1],linewidth=1)
        end
    end
    if cb
        #cl2       = nc-mc > 0 ? nc-mc : 2
        color_bar = heatmap(rand(2,2),clims=(mc,nc+1),framestyle=:none,c=pal,cbar=true,lims=(-1,0))  
        return plot(plt,color_bar,layout=l)
    else
        return plt
    end
end;
