struct SeenEdge{I<:Integer}
    nodes::SVector{2,I}
    index::I #stores the index of the node generated when bisecting the edge 
end
SeenEdge(v::SVector{2,I},i::J) where {I<:Integer,J<:Integer} = SeenEdge(v,convert(I,i))
@inline Base.hash(sedge::SeenEdge) = hash(sedge.nodes)
@inline Base.isequal(v::SVector{2,I},se::SeenEdge{I}) where I<:Integer = isequal(v,se.nodes) || isequal(reverse(v),se.nodes)
@inline Base.isequal(se::SeenEdge{I},v::SVector{2,I}) where I<:Integer = isequal(v,se)

struct RefAux{I<:Integer,P<:Integer}
    degs::MVector{6,P}
    dots::MVector{6,I}
    seen::FESet{SeenEdge{I}}
end
RefAux{I,P}() where {I<:Integer,P<:Integer} = RefAux(MVector{6,P}(zeros(6)),MVector{6,I}(zeros(6)),FESet{SeenEdge{I}}())

function index(seen::FESet{SeenEdge{I}},v::SVector{2,I}) where I<:Integer
    isin,(_,k) = gettoken(seen,v)
    if isin
        return gettokenvalue(seen,k).index
    else
        isin,(_,k) = gettoken(seen,reverse(v))
        if isin
            return gettokenvalue(seen,k).index
        else
            return 0
        end
    end
end
 """

    mark_triangles(estim,mesh::HPMesh;estim_params = nothing)

Marks triangles with vertices `vert` for which `estim(vert,estim_param)` returns `true`. It only performs the marking of the triangles (in `mesh.rT`) and of its edges, in `mesh.rE`. In order to mark further for preserving the conformity of the mesh it is necessary to run `h_comformity!(mesh)`.
"""
function mark_triangles!(estim,mesh::MeshHP;estim_params...)
    (;points,triangles,edges) = mesh
    tri  = MMatrix{2,3}(zeros(2,3))
    for t in triangles
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
        for t in triangles 
            num_marked = count(ismarked,edges[get_edges(t)])     
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


function refine_red!(mesh::MeshHP{F,I,P,Bool},refaux::RefAux{I,P},i) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edges,triangles) = mesh
    (;degs,dots,seen) = refaux
    for t in filter(isred,triangles)
        dots[1:3] .= nodes(t)
        t_edges    = get_edges(t) 
        degs[1:3] .= [degree(edges[e]) for e in t_edges]
        degs[4:6] .= degs[1:3]-degs[[3,1,2]] |> x->max.(abs.(x),1)
        for j in eachindex(t_edges)
            edge = t_edges[j]
            k    = index(seen,edge)
            if k>0
                dots[j+3] = k
                remove!(edges,edges[edge])
            else
                points[:,i] = sum(points[:,edge],dims=2)/2.
                dots[j+3]   = i
                add!(seen,SeenEdge(edge,i))
                oldedge = edges[edge]
                m = marker(oldedge)
                add!(edges,EdgeHP{I,P}([edge[1],i],degs[j],m,false))
                add!(edges,EdgeHP{I,P}([i,edge[2]],degs[j],m,false))
                i += 1
            end
        end 
        add!(edges,EdgeHP{I,P}([dots[6],dots[4]],degs[4],0,false))
        add!(edges,EdgeHP{I,P}([dots[4],dots[5]],degs[5],0,false))    
        add!(edges,EdgeHP{I,P}([dots[5],dots[6]],degs[6],0,false))
        add!(triangles,TriangleHP{I,P}(dots[[1,4,6]],0))
        add!(triangles,TriangleHP{I,P}(dots[[4,2,5]],0))
        add!(triangles,TriangleHP{I,P}(dots[[6,5,3]],0))
        add!(triangles,TriangleHP{I,P}(dots[[5,6,4]],0))
        remove!(triangles,t)
    end
    return i
end

function refine_blue!(mesh::MeshHP{F,I,P,Bool},refaux::RefAux{I,P},i) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edges,triangles) = mesh
    (;degs,dots,seen) = refaux
    for t in filter(isblue,triangles)
        dots[1:3] = nodes(t)
        t_edges   = get_edges(t)
        degs[1:3] .= [degree(edges[e]) for e in get_edges(t)]
        degs[4]    = max(1,abs(degs[1]-degs[2]),abs(degs[1]-degs[3]))
        if ismarked(edges[t_edges[2]])
            degs[5] = max(1,abs(degs[1]-degs[2]),abs(degs[2]-degs[4]))
            for j in 1:2
                edge = t_edges[j]
                k    = index(seen,edge)
                if k>0
                    dots[j+3] = k
                    remove!(edges,edges[edge])
                else
                    points[:,i] = sum(points[:,edge],dims=2)/2.
                    dots[j+3]   = i
                    add!(seen,SeenEdge(edge,i))
                    oldedge = edges[edge]
                    m = marker(oldedge)
                    add!(edges,EdgeHP{I,P}([edge[1],i],degs[j],m,false))
                    add!(edges,EdgeHP{I,P}([i,edge[2]],degs[j],m,false))
                    i += 1
                end
            end
            add!(edges,EdgeHP{I,P}([dots[3],dots[4]],degs[4],0,false))
            add!(edges,EdgeHP{I,P}([dots[5],dots[4]],degs[5],0,false))    
            add!(triangles,TriangleHP{I,P}(dots[[1,4,3]],0))
            add!(triangles,TriangleHP{I,P}(dots[[4,2,5]],0))
            add!(triangles,TriangleHP{I,P}(dots[[4,5,3]],0))
            remove!(triangles,t)
        elseif ismarked(edges[t_edges[3]])
            degs[5] = max(1,abs(degs[1]-degs[3]),abs(degs[3]-degs[4]))
            for j in 0:1
                edge = t_edges[1+2j]
                k    = index(seen,edge)
                if k>0
                    dots[j+4] = k
                    remove!(edges,edges[edge])
                else
                    points[:,i] = sum(points[:,edge],dims=2)/2.
                    dots[j+4]   = i
                    add!(seen,SeenEdge(edge,i))
                    oldedge = edges[edge]
                    m = marker(oldedge)
                    add!(edges,EdgeHP{I,P}([edge[1],i],degs[1+2j],m,false))
                    add!(edges,EdgeHP{I,P}([i,edge[2]],degs[1+2j],m,false))
                    i += 1
                end
            end
            add!(edges,EdgeHP{I,P}([dots[3],dots[4]],degs[4],0,false))
            add!(edges,EdgeHP{I,P}([dots[4],dots[5]],degs[5],0,false))    
            add!(triangles,TriangleHP{I,P}(dots[[1,4,5]],0))
            add!(triangles,TriangleHP{I,P}(dots[[4,2,3]],0))
            add!(triangles,TriangleHP{I,P}(dots[[4,3,5]],0))
            remove!(triangles,t)
        end
    end
    return i
end
function refine_green!(mesh::MeshHP{F,I,P,Bool},refaux::RefAux{I,P},i) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edges,triangles) = mesh
    (;degs,dots,seen) = refaux
    for t in filter(isgreen,triangles)
        dots[1:3] = nodes(t)
        edge = longestedge(t)
        degs[1:3] .= [degree(edges[e]) for e in get_edges(t)]
        degs[4]   = max(1,abs(degs[1]-degs[2]),abs(degs[1]-degs[3]))
        k = index(seen,edge)
        if k>0
            dots[4] = k
            remove!(edges,edges[edge])
        else
            points[:,i] = sum(points[:,edge],dims=2)/2.
            dots[4]     = i
            add!(seen,SeenEdge(edge,i))
            oldedge = edges[edge] 
            m = marker(oldedge)
            add!(edges,EdgeHP{I,P}([dots[1],dots[4]],degs[1],0,false))
            add!(edges,EdgeHP{I,P}([dots[4],dots[2]],degs[1],0,false))
            i  += 1
        end
        add!(edges,EdgeHP{I,P}([dots[3],dots[4]],degs[4],0,false))
        add!(triangles,TriangleHP{I,P}(dots[[1,4,3]],0))
        add!(triangles,TriangleHP{I,P}(dots[[4,2,3]],0))
        remove!(triangles,t)
    end
    return i
end

    
function refine!(mesh::MeshHP{F,I,P,Bool}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edges,triangles) = mesh
    i = size(points,2) + 1
    
    append!(points,zeros(F,2,count(ismarked,edges)))
    refaux  = RefAux{I,P}()
    i = refine_red!(mesh,refaux,i)
    i = refine_blue!(mesh,refaux,i)
    i = refine_green!(mesh,refaux,i)
    for edge in edges
        if ismarked(edge)
            remove!(edges,edge)
        end
    end
end
