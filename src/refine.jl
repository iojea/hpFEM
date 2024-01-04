# struct SeenEdge{I<:Integer}
#     nodes::SVector{2,I}
#     index::I #stores the index of the node generated when bisecting the edge 
# end
# SeenEdge(v::SVector{2,I},i::J) where {I<:Integer,J<:Integer} = SeenEdge(v,convert(I,i))
# @inline Base.hash(sedge::SeenEdge) = hash(sedge.nodes)
# @inline Base.isequal(v::SVector{2,I},se::SeenEdge{I}) where I<:Integer = isequal(v,se.nodes) || isequal(reverse(v),se.nodes)
# @inline Base.isequal(se::SeenEdge{I},v::SVector{2,I}) where I<:Integer = isequal(v,se)

struct RefAux{I<:Integer,P<:Integer}
    degs::MVector{6,P}
    dots::MVector{6,I}
    seen::DictHP{2,I,I}
end
RefAux{I,P}() where {I<:Integer,P<:Integer} = RefAux(MVector{6,P}(zeros(6)),MVector{6,I}(zeros(6)),DictHP{2,I,I}())
function RefAux(mesh::MeshHP{F,I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    seen = fill(I(0),filter(ismarked,mesh.edges))
    RefAux(MVector{6,P}(zeros(6)),MVector{6,I}(zeros(6)),seen)
end

# function index(seen::FESet{SeenEdge{I}},v::SVector{2,I}) where I<:Integer
#     isin,(_,k) = gettoken(seen,v)
#     if isin
#         return gettokenvalue(seen,k).index
#     else
#         isin,(_,k) = gettoken(seen,reverse(v))
#         if isin
#             return gettokenvalue(seen,k).index
#         else
#             return 0
#         end
#     end
# end
 """

    mark_triangles(estim,mesh::HPMesh;estim_params = nothing)

Marks triangles with vertices `vert` for which `estim(vert,estim_param)` returns `true`. It only performs the marking of the triangles (in `mesh.rT`) and of its edges, in `mesh.rE`. In order to mark further for preserving the conformity of the mesh it is necessary to run `h_comformity!(mesh)`.
"""
function mark_triangles!(estim,mesh::MeshHP;estim_params...)
    (;points,triangles,edges) = mesh
    tri  = MMatrix{2,3}(zeros(2,3))
    for t in pairs(triangles)
        tri .= points[:,nodes(t)]            
        #if estim(@views points[:,SVector(triangle(t))];estim_params...)
        if estim(tri;estim_params...) 
            mark!(t)
            for e in get_edges(t)
                mark!(edges[e])
            end
        end  
    end
end

function estim_distance_origin(vert;h=0.2,μ=1)
    ℓ = minimum(norm(vert[:,i]-vert[:,j]) for (i,j) in ((1,2),(1,3),(2,3)))
    d = maximum(norm(v) for v in eachcol(vert))
    if intriangle(zeros(2),vert)≥0
        return ℓ>h^(1/μ)
    else
        return ℓ>1.5h*d^(1-μ) 
    end    
end

function intriangle(p::T,a::V,b::V,c::V) where {T<:AbstractArray,V<:AbstractArray}
    if abs(orient(a,b,p)+orient(b,c,p)+orient(c,a,p))==3
        return 1
        elseif orient(a,b,p)*orient(b,c,p)*orient(c,a,p) == 0
        return 0
    else
        return -1
    end
end
intriangle(p::T,vert::M) where {T<:AbstractArray,M<:AbstractArray} = intriangle(p,eachcol(vert)...)



"""

    h_conformity!(mesh::HPMesh)

Performs the marking of previously un-marked triangles (and corresponding edges) in order to avoid hanging nodes. This function must be run after the marking of the unsuitable triangles, by `mark_triangles`.
"""

function h_conformity!(mesh::MeshHP)
    (;points,triangles,edges) = mesh
    still = true
    while still
        still = false
        for t in pairs(triangles)
            num_marked = count(ismarked,getindices(edges,get_edges(t)))     
            if num_marked>0
                long_edge = edges[longestedge(t)]
                if !ismarked(long_edge)
                    mark!(long_edge)
                    still = true
                    num_marked += 1
                end
                mark!(t,num_marked)
            end
        end
    end
end


function refine_red!(mesh::MeshHP{F,I,P},refaux::RefAux{I,P},i) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edges,triangles) = mesh
    (;degs,dots,seen) = refaux
    for t in filter(isred,pairs(triangles))
        dots[1:3] .= nodes(t)
        t_edges    = get_edges(t) 
        degs[1:3] .= [degree(edges[e]) for e in t_edges]
        degs[4:6] .= degs[1:3]-degs[[3,1,2]] |> x->max.(abs.(x),1)
        for j in eachindex(t_edges)
            edge = t_edges[j]
            k    = seen[edge]
            if k>0
                dots[j+3] = k
                remove!(edges,edge)
            else
                points[:,i] = sum(points[:,edge],dims=2)/2.
                dots[j+3]   = i
                add!(seen,edge,i)
                oldedge = edges[edge]
                m = marker(oldedge)
                add!(edges,EdgeHP{I}(edge[1],i),EdgeProperties{P,Bool}(degs[j],m,false))
                add!(edges,EdgeHP{I}(i,edge[2]),EdgeProperties{P,Bool}(degs[j],m,false))
                i += 1
            end
        end 
        add!(edges,EdgeHP{I}(dots[6],dots[4]),EdgeProperties{P,Bool}(degs[4],0,false))
        add!(edges,EdgeHP{I}(dots[4],dots[5]),EdgeProperties{P,Bool}(degs[5],0,false))    
        add!(edges,EdgeHP{I}(dots[5],dots[6]),EdgeProperties{P,Bool}(degs[6],0,false))
        add!(triangles,TriangleHP{I}(dots[[1,4,6]]),TriangleProperties{P}(0))
        add!(triangles,TriangleHP{I}(dots[[4,2,5]]),TriangleProperties{P}(0))
        add!(triangles,TriangleHP{I}(dots[[6,5,3]]),TriangleProperties{P}(0))
        add!(triangles,TriangleHP{I}(dots[[5,6,4]]),TriangleProperties{P}(0))
        remove!(triangles,nodes(t))
    end
    return i
end

function refine_blue!(mesh::MeshHP{F,I,P},refaux::RefAux{I,P},i) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edges,triangles) = mesh
    (;degs,dots,seen) = refaux
    for t in filter(isblue,pairs(triangles))
        dots[1:3] = nodes(t)
        t_edges   = get_edges(t)
        degs[1:3] .= [degree(edges[e]) for e in t_edges]
        degs[4]    = max(1,abs(degs[1]-degs[2]),abs(degs[1]-degs[3]))
        if ismarked(edges[t_edges[2]])
            degs[5] = max(1,abs(degs[1]-degs[2]),abs(degs[2]-degs[4]))
            for j in 1:2
                edge = t_edges[j]
                k    = seen[edge]
                if k>0
                    dots[j+3] = k
                    remove!(edges,edge)
                else
                    points[:,i] = sum(points[:,edge],dims=2)/2.
                    dots[j+3]   = i
                    add!(seen,edge,i)
                    oldedge = edges[edge]
                    m = marker(oldedge)
                    add!(edges,EdgeHP{I}(edge[1],i),EdgeProperties{P,Bool}(degs[j],m,false))
                    add!(edges,EdgeHP{I}(i,edge[2]),EdgeProperties{P,Bool}(degs[j],m,false))
                    i += 1
                end
            end
            add!(edges,EdgeHP{I}(dots[3],dots[4]),EdgeProperties{P,Bool}(degs[4],0,false))
            add!(edges,EdgeHP{I}(dots[5],dots[4]),EdgeProperties{P,Bool}(degs[5],0,false))    
            add!(triangles,TriangleHP{I}(dots[[1,4,3]]),TriangleProperties{P}(0))
            add!(triangles,TriangleHP{I}(dots[[4,2,5]]),TriangleProperties{P}(0))
            add!(triangles,TriangleHP{I}(dots[[4,5,3]]),TriangleProperties{P}(0))
            remove!(triangles,nodes(t))
        elseif ismarked(edges[t_edges[3]])
            degs[5] = max(1,abs(degs[1]-degs[3]),abs(degs[3]-degs[4]))
            for j in 0:1
                edge = t_edges[1+2j]
                k    = seen[edge]
                if k>0
                    dots[j+4] = k
                    remove!(edges,edge)
                else
                    points[:,i] = sum(points[:,edge],dims=2)/2.
                    dots[j+4]   = i
                    add!(seen,edge,i)
                    oldedge = edges[edge]
                    m = marker(oldedge)
                    add!(edges,EdgeHP{I}(edge[1],i),EdgeProperties{P,Bool}(degs[1+2j],m,false))
                    add!(edges,EdgeHP{I}(i,edge[2]),EdgeProperties{P,Bool}(degs[1+2j],m,false))
                    i += 1
                end
            end
            add!(edges,EdgeHP{I}(dots[3],dots[4]),EdgeProperties{P,Bool}(degs[4],0,false))
            add!(edges,EdgeHP{I}(dots[4],dots[5]),EdgeProperties{P,Bool}(degs[5],0,false))    
            add!(triangles,TriangleHP{I}(dots[[1,4,5]]),TriangleProperties{P}(0))
            add!(triangles,TriangleHP{I}(dots[[4,2,3]]),TriangleProperties{P}(0))
            add!(triangles,TriangleHP{I}(dots[[4,3,5]]),TriangleProperties{P}(0))
            remove!(triangles,nodes(t))
        end
    end
    return i
end
function refine_green!(mesh::MeshHP{F,I,P},refaux::RefAux{I,P},i) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edges,triangles) = mesh
    (;degs,dots,seen) = refaux
    for t in filter(isgreen,pairs(triangles))
        dots[1:3] = nodes(t)
        edge = longestedge(t)
        degs[1:3] .= [degree(edges[e]) for e in get_edges(t)]
        degs[4]   = max(1,abs(degs[1]-degs[2]),abs(degs[1]-degs[3]))
        k = seen[edge]
        if k>0
            dots[4] = k
            remove!(edges,edge)
        else
            points[:,i] = sum(points[:,edge],dims=2)/2.
            dots[4]     = i
            add!(seen,edge,i)
            oldedge = edges[edge] 
            m = marker(oldedge)
            add!(edges,EdgeHP{I}(dots[1],dots[4]),EdgeProperties{P,Bool}(degs[1],0,false))
            add!(edges,EdgeHP{I}(dots[4],dots[2]),EdgeProperties{P,Bool}(degs[1],0,false))
            i  += 1
        end
        add!(edges,EdgeHP{I}(dots[3],dots[4]),EdgeProperties{P,Bool}(degs[4],0,false))
        add!(triangles,TriangleHP{I}(dots[[1,4,3]]),TriangleProperties{P}(0))
        add!(triangles,TriangleHP{I}(dots[[4,2,3]]),TriangleProperties{P}(0))
        remove!(triangles,nodes(t))
    end
    return i
end

    
function refine!(mesh::MeshHP{F,I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edges,triangles) = mesh
    i = size(points,2) + 1
    println(size(points))
    append!(points,zeros(F,2,count(ismarked,edges)))
    println(size(points))
    refaux  = RefAux(mesh)
    println(refaux.seen)
    i = refine_red!(mesh,refaux,i)
    println(refaux.seen)
    i = refine_blue!(mesh,refaux,i)
    println(refaux.seen)
    i = refine_green!(mesh,refaux,i)
    println(refaux.seen)
    filter!(!ismarked,mesh.edges)
    # for edge in pairs(edges)
    #     if ismarked(edge)
    #         unset!(edges,nodes(edge))
    #     end
    # end
end
