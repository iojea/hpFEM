struct RefAux{I<:Integer,P<:Integer}
    i::Ref{I}
    degs::MVector{6,P}
    dots::MVector{6,I}
    seen::Dictionary{EdgeHP{I},I}
end
RefAux{I,P}() where {I<:Integer,P<:Integer} = RefAux(MVector{6,P}(zeros(6)),MVector{6,I}(zeros(6)),Dictionary{EdgeHP{I},I}())
function RefAux(i,mesh::MeshHP{F,I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    seen = fill(I(0),filter(ismarked,mesh.edgelist))
    RefAux(Ref(i),MVector{6,P}(zeros(6)),MVector{6,I}(zeros(6)),seen)
end

 """

    mark!(estim,mesh::HPMesh;estim_params...)

Marks triangles with vertices `vert` for which `estim(vert,estim_param)` returns `true`. It only performs the marking of the triangles (in `mesh.rT`) and of its edgelist, in `mesh.rE`. In order to mark further for preserving the conformity of the mesh it is necessary to run `h_comformity!(mesh)`.
"""
function mark!(estim,mesh::MeshHP;estim_params...)
    (;points,trilist,edgelist) = mesh
    tri  = MMatrix{2,3}(zeros(2,3))
    for t in triangles(trilist)
        tri .= points[:,t]            
        if estim(tri;estim_params...) 
            for e in edges(t)
                mark!(edgelist[e])
            end
        end  
    end
    _h_conformity!(mesh)
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

Performs the marking of previously un-marked trilist (and corresponding edges) in order to avoid hanging nodes.
"""

function _h_conformity!(mesh::MeshHP)
    (;trilist,edgelist) = mesh
    still = true
    while still
        still = false
        for t in triangles(trilist)
            num_marked = count(ismarked,getindices(edgelist,edges(t)))     
            if num_marked>0
                long_edge = edgelist[longestedge(t)]
                if !ismarked(long_edge)
                    mark!(long_edge)
                    still = true
                    num_marked += 1
                end
                mark!(trilist[t],num_marked)
            end
        end
    end
end


function refine_red!(t::TriangleHP{I},mesh::MeshHP{F,I,P},refaux::RefAux{I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edgelist,trilist) = mesh
    (;i,degs,dots,seen) = refaux
    dots[1:3] .= t
    t_edges    = edges(t) 
    degs[1:3] .= [degree(edgelist[e]) for e in t_edges]
    degs[4:6] .= degs[1:3]-degs[[3,1,2]] |> x->max.(abs.(x),1)
    for j in eachindex(t_edges)
        edge = t_edges[j]
        k    = seen[edge]
        if k>0
            dots[j+3] = k
        else
            points[:,i[]] = sum(points[:,edge],dims=2)/2.
            dots[j+3]   = i[]
            set!(seen,edge,i[])
            m = marker(edgelist[edge])
            set!(edgelist,EdgeHP{I}(edge[1],i[]),EdgeProperties{P,Bool}(degs[j],m,false))
            set!(edgelist,EdgeHP{I}(i[],edge[2]),EdgeProperties{P,Bool}(degs[j],m,false))
            i[] += 1
        end
    end 
    set!(edgelist,EdgeHP{I}(dots[6],dots[4]),EdgeProperties{P,Bool}(degs[4],0,false))
    set!(edgelist,EdgeHP{I}(dots[4],dots[5]),EdgeProperties{P,Bool}(degs[5],0,false))    
    set!(edgelist,EdgeHP{I}(dots[5],dots[6]),EdgeProperties{P,Bool}(degs[6],0,false))
    set!(trilist,TriangleHP{I}(dots[[1,4,6]]),TriangleProperties{P}())
    set!(trilist,TriangleHP{I}(dots[[4,2,5]]),TriangleProperties{P}())
    set!(trilist,TriangleHP{I}(dots[[6,5,3]]),TriangleProperties{P}())
    set!(trilist,TriangleHP{I}(dots[[5,6,4]]),TriangleProperties{P}())
end

function refine_blue!(t::TriangleHP{I},mesh::MeshHP{F,I,P},refaux::RefAux{I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edgelist,trilist) = mesh
    (;i,degs,dots,seen) = refaux
    dots[1:3] .= t
    t_edges    = edges(t)
    degs[1:3] .= [degree(edgelist[e]) for e in t_edges]
    degs[4]    = max(1,abs(degs[1]-degs[2]),abs(degs[1]-degs[3]))
    if ismarked(edgelist[t_edges[2]])
        degs[5] = max(1,abs(degs[1]-degs[2]),abs(degs[2]-degs[4]))
        for j in 1:2
            edge = t_edges[j]
            k    = seen[edge]
            if k>0
                dots[j+3] = k
            else
                points[:,i[]] = sum(points[:,edge],dims=2)/2.
                dots[j+3]   = i[]
                set!(seen,edge,i[])
                m = marker(edgelist[edge])
                set!(edgelist,EdgeHP{I}(edge[1],i[]),EdgeProperties{P,Bool}(degs[j],m,false))
                set!(edgelist,EdgeHP{I}(i[],edge[2]),EdgeProperties{P,Bool}(degs[j],m,false))
                i[] += 1
            end
        end
        set!(edgelist,EdgeHP{I}(dots[3],dots[4]),EdgeProperties{P,Bool}(degs[4],0,false))
        set!(edgelist,EdgeHP{I}(dots[5],dots[4]),EdgeProperties{P,Bool}(degs[5],0,false))    
        set!(trilist,TriangleHP{I}(dots[[1,4,3]]),TriangleProperties{P}())
        set!(trilist,TriangleHP{I}(dots[[4,2,5]]),TriangleProperties{P}())
        set!(trilist,TriangleHP{I}(dots[[4,5,3]]),TriangleProperties{P}())
    elseif ismarked(edgelist[t_edges[3]])
        degs[5] = max(1,abs(degs[1]-degs[3]),abs(degs[3]-degs[4]))
        for j in 0:1
            edge = t_edges[1+2j]
            k    = seen[edge]
            if k>0
                dots[j+4] = k
            else
                points[:,i[]] = sum(points[:,edge],dims=2)/2.
                dots[j+4]   = i[]
                set!(seen,edge,i[])
                m = marker(edgelist[edge])
                set!(edgelist,EdgeHP{I}(edge[1],i[]),EdgeProperties{P,Bool}(degs[1+2j],m,false))
                set!(edgelist,EdgeHP{I}(i[],edge[2]),EdgeProperties{P,Bool}(degs[1+2j],m,false))
                i[] += 1
            end
        end
        set!(edgelist,EdgeHP{I}(dots[3],dots[4]),EdgeProperties{P,Bool}(degs[4],0,false))
        set!(edgelist,EdgeHP{I}(dots[4],dots[5]),EdgeProperties{P,Bool}(degs[5],0,false))    
        set!(trilist,TriangleHP{I}(dots[[1,4,5]]),TriangleProperties{P}())
        set!(trilist,TriangleHP{I}(dots[[4,2,3]]),TriangleProperties{P}())
        set!(trilist,TriangleHP{I}(dots[[4,3,5]]),TriangleProperties{P}())
    end
end


function refine_green!(t::TriangleHP{I},mesh::MeshHP{F,I,P},refaux::RefAux{I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edgelist,trilist) = mesh
    (;i,degs,dots,seen) = refaux
    dots[1:3] .= t
    edge = longestedge(t)
    degs[1:3] .= [degree(edgelist[e]) for e in edges(t)]
    degs[4]   = max(1,abs(degs[1]-degs[2]),abs(degs[1]-degs[3]))
    k = seen[edge]
    if k>0
        dots[4] = k
    else
        points[:,i[]] = sum(points[:,edge],dims=2)/2.
        dots[4]     = i[]
        set!(seen,edge,i[])
        oldedge = edgelist[edge] 
        m = marker(oldedge)
        set!(edgelist,EdgeHP{I}(dots[1],dots[4]),EdgeProperties{P,Bool}(degs[1],0,false))
        set!(edgelist,EdgeHP{I}(dots[4],dots[2]),EdgeProperties{P,Bool}(degs[1],0,false))
        i[] += 1
    end
    set!(edgelist,EdgeHP{I}(dots[3],dots[4]),EdgeProperties{P,Bool}(degs[4],0,false))
    set!(trilist,TriangleHP{I}(dots[[1,4,3]]),TriangleProperties{P}())
    set!(trilist,TriangleHP{I}(dots[[4,2,3]]),TriangleProperties{P}())
end

    
function refine!(mesh::MeshHP{F,I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edgelist,trilist) = mesh
    i = I(size(points,2) + 1)
    n_edgelist = count(ismarked,edgelist)
    n_tris  = 4count(isred,trilist)+3count(isblue,trilist)+2count(isgreen,trilist)
    sizehint!(trilist,length(trilist)+n_tris)
    sizehint!(edgelist,length(trilist)+n_edgelist)
    append!(points,zeros(F,2,n_edgelist))
    refaux  = RefAux(i,mesh)
    for t in triangles(trilist)
        if isred(trilist[t])
            refine_red!(t,mesh,refaux)
        elseif isblue(trilist[t])
            refine_blue!(t,mesh,refaux)
        elseif isgreen(trilist[t])
            refine_green!(t,mesh,refaux)
        end
    end
    filter!(!ismarked,mesh.trilist)
    filter!(!ismarked,mesh.edgelist)
end



function check_p_conformity(pt::Vector{T}) where T<:Integer
    sum(pt) - maximum(pt) ≥ maximum(pt) 
end

function p_conformity!(mesh::MeshHP{F,I,P},t::TriangleHP{I},d) where {F,I,P}
    (;points,trilist,edgelist) = mesh
    p,eds = pedges(t)
    out = false
    if check_p_conformity(p)
        out = true
    else
        if d>0
            setdegree!(edgelist[eds[1]],p[1] + 1)
            t₁ = neighbor(mesh,t,eds[1])
            if p_conformity!(mesh,t₁,d-1)
                out = true
            else
                setdegree!(edgelist[eds[1]],p[1])
                setdegree!(edgelist[eds[2]],p[2]+1)
                t₂ = neighbor(mesh,t,eds[2])
                if p_conformity!(mesh,t₂,d-1)
                    out = true
                else
                    setdegree!(edgelist[eds[2]],p[2])
                end
            end
        end
    end
    out
end
function p_conformity!(mesh::MeshHP{F,I,P}) where {F,I,P}
    (;trilist) = mesh
    for t in triangles(trilist)
        println(t)
        p_conformity!(mesh,t,10)
    end
end

function neighbor(mesh::MeshHP{F,I,P},t::TriangleHP{I},e::EdgeHP{I}) where {F,I,P}
    for t′ in triangles(mesh)
        if t′!=t && (e in edges(t′))
            return t′
        end
    end
    return nothing
end
