"""
    normals(mesh::MeshHP)  

returns a dictionary with keys given by `edges(mesh)` and values of type `SVector{2,Float64}` with the normal vectors to each edge.
"""
function normals(mesh::MeshHP)
    (; points, edgelist) = mesh
    dict = similar(edgelist, SVector{2,Float64})
    for e in edges(edgelist)
        pts = view(points, :, e)
        normal = [pts[2, 1] - pts[2, 2], pts[1, 2] - pts[1, 1]]
        set!(dict, e, normal / norm(normal))
    end
    dict
end


function find_neighbors(mesh::MeshHP{F,I,P}) where {F,I,P}
    (; edgelist, trilist) = mesh
    neighs = map(copy, fill(Vector{TriangleHP{I}}(), edgelist))
    for t in keys(trilist)
        vecs = getindices(neighs, edges(t))
        (x -> push!(x, t)).(vecs)
    end
    neighs
end


function jumps(mesh::MeshHP, bd::BasisDict, dof, sol)
    (; edgelist) = mesh
    nrms = normals(mesh)
    jmp = similar(nrms)
    for e in keys(filter(isinterior, edgelist))
        pₑ = degree(edgelist[e])
        v = adjacent(edgelist[e])
        t1 = (e[1], e[2], v[1])
        t2 = (e[1], e[2], v[2])
        dof1 = dof[t1]
        dof2 = dof[t2]
        p1, _ = pnodes(t1, mesh)
        p2, _ = pnodes(t2, mesh)
        ∇b1 = bd[p1].∇b
        ∇b2 = bd[p2].∇b
        sch = grundmann_moeller(Float64, Val(1), 4 * pₑ + 1)
        lin(t) = t * e[2] + (1 - t) * e[1]
        ∂b1∂ηₑ = t -> sum(sol[dof1[k]] * ∇b1[k](lin(t)) ⋅ nrms[e] for k in eachindex(dof1))
        ∂b2∂ηₑ = t -> sum(sol[dof2[k]] * ∇b2[k](lin(t)) ⋅ nrms[e] for k in eachindex(dof2))
        jmp[e] = integrate(t -> (∂b1∂ηₑ - ∂b2∂ηₑ)^2, sch, [0, 1])
    end
    jmp
end


h(vert) = maximum(map(norm,(vert[:,1]-vert[:,2],vert[:,2]-vert[:,3],vert[:,3]-vert[:,1])))
h(mesh::MeshHP,t) = h(mesh.trilist[:,t])


D(vert) = maximum(norm.(eachcol(vert)))
D(mesh::MeshHP,t) = D(mesh.trilist[:,t])


#function η(mesh::MeshHP)
#    
#end