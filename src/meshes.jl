"""
    MeshHP{F<:AbstractFloat,I<:Integer,P<:Integer}()

A mesh for `HP` finite element methods. Its fields are: 
    + `points::Matrix{F}`: a matrix of points
    + `triangles::FESet{TriangleHP{I,P,B}}`: set of triangles
    + `edges::FESet{EdgeHP{I,P,B}}`: set of edges

    MeshHP(tri::TriangulationIO)
builds an `MeshHP` from a `Triangulate.TriangulatioIO` struct.  
"""
struct MeshHP{F<:AbstractFloat,I<:Integer,P<:Integer,Bool}
    points::ElasticMatrix{F,Vector{F}}
    triangles::FESet{TriangleHP{I,P}}
    edges::FESet{EdgeHP{I,P,Bool}}    
end
MeshHP(mat::Matrix,tris,edgs) = MeshHP(ElasticMatrix(mat),tris,edgs)
Base.copy(mesh::MeshHP) = MeshHP(deepcopy(mesh.points),deepcopy(mesh.triangles),deepcopy(mesh.edges))

@inline get_edges(t::TriangleHP,mesh::MeshHP) = [mesh.edges[e] for e in get_edges(t)]

function MeshHP(tri::TriangulateIO)
    (;pointlist,trianglelist,edgelist,edgemarkerlist) = tri
    F = eltype(pointlist)
    I = eltype(trianglelist)
    P = UInt8
    triangles = FESet([triangle(v,pointlist) for v in eachcol(trianglelist)])
    edges = FESet([EdgeHP{I,P,Bool}(e,1,edgemarkerlist[i],false) for (i,e) in enumerate(eachcol(edgelist))] )
    MeshHP(pointlist,triangles,edges)
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
    return MeshHP(tri)
end;
circular_mesh(r,h) = circular_mesh(zeros(2),r,h)
circular_mesh(h)   = circular_mesh(1,h)


function correct_boundary_circular(mesh::MeshHP)
    (;points,edges) = mesh
    for e in edges
        if marker(e)==1
            i,j = nodes(e)
            points[:,i] .= points[:,i]/norm(points[:,i])
            points[:,j] .= points[:,j]/norm(points[:,j])
        end
    end
end

function circular_mesh_graded_to_center(h,μ;maxiter=4,rec=false)
    k = 0
    mesh = circular_mesh(h)
    flag = true
    mshs = [copy(mesh)]
    while k≤maxiter && flag
        mark_triangles!(estim_distance_origin,mesh,h=h,μ=μ)
        rec ? push!(mshs,copy(mesh)) : nothing
        h_conformity!(mesh)
        rec ? push!(mshs,copy(mesh)) : nothing
        if count(ismarked,mesh.triangles) > 0
            refine!(mesh)
            correct_boundary_circular(mesh)
            k += 1
            println("Iteration:",k)
            rec ? push!(mshs,copy(mesh)) : nothing
        else
            flag = false
        end
    end
    return rec ? mshs : mesh
end



function psortperm(v)
    i    = argmin(v)
    if v[i] ≤ v[mod1(i+1,3)] ≤ v[mod1(i+2,3)]
        ind = mod1.(i:i+2,3)
    else
        ind = mod1.(i:-1:i-2,3)
    end
    return ind
end
