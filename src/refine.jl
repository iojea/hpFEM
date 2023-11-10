include("meshes.jl")

#"""
#    mutable struct RefInd{T}
#
#a structure for storing the indices necessary for tracking the construction of the refined mesh
#"""
struct RefInd{U<:Unsigned,N,N1,N2,V0<:MVector{3}{U},V<:MVector{N}{U},V1<:MVector{N1}{U},V2<:MVector{N2}{U}} 
    PTE::V0 #vector of indices for next Point, Triangle and Edge. 
    PE::V #Vector storing the index of the node created by bisecting an edge 
    EE::V #Vector storing the index of the first edge that is a child of the current edge
    nodes::V1 #Stores the indices of the new nodes (at most 6)
    edges::V2 #Stores the indices of the new edges (at most 9)
end


cross2(u,v) = u[1]*v[2]-u[2]*v[1]
function in_triangle(point,vert)
    i1 = argmin(norm(v) for v in eachcol(vert))
    i2 = mod1(i1+1,3)
    i3 = mod1(i1+2,3)
    a = (cross2(point,vert[:,i3])-cross2(vert[:,i1],vert[:,i3]))/cross2(vert[:,i2],vert[:,i3])
    b = (-cross2(point,vert[:,i2])+cross2(vert[:,i1],vert[:,i2]))/cross2(vert[:,i2],vert[:,i3])
    a≥0 && b≥0 && a+b≤1
end
    
"""

    compute_leT!(mesh::HPMesh)

Computes the largest edge of each triangle, and stores it in `mesh.leT`. This functions initializes the mesh for refining, since the longest edge is the "marked" edge. 
"""
function compute_leT!(mesh::HPMesh)
    (;P,ET,leT) = mesh
    for k in eachindex(eachcol(mesh.T))
        edges  = view(mesh.E,:,ET[:,k])
        leT[k] = UInt8(argmax(norm(P[:,edges[1,j]]-P[:,edges[2,j]]) for j in 1:3))
    end
end


"""

    mark_triangles(estima,mesh::HPMesh;estim_params = nothing)

Marks triangles with vertices `vert` for which `estim(vert,estim_param)` returns `true`. It only performs the marking of the triangles (in `mesh.rT`) and of its edges, in `mesh.rE`. In order to mark further for preserving the conformity of the mesh it is necessary to run `h_comformity!(mesh)`.
"""
function mark_triangles!(estim,mesh::HPMesh;estim_params...)
    (;P,T,ET,rE,rT) = mesh
    #rE = view(mesh.rE,:)
    #rT = view(mesh.rT,:)
    #ET = view(mesh.ET,:,:)
    for k in eachindex(rT)
        if estim(@views P[:,T[:,k]];estim_params...)
            # rT[k]   = 4
            rE[ET[:,k]] .= true
        end  
    end
end

"""

    codify_ref(a,b)

codifies two integers `1≤a≠b≤3` as `c=4k+a`, where `k+a = b`. It is the code for storing two edges to be refined as a single integer.  
"""
codify_ref(a,b) = 4((3+b-a)%3)+a 

"""

    decodify_ref(c)

decodes the `c=4k+a` as two different integers `1≤a,b≤3` that correspond two the indices of two edges to be refined. 
"""
function decodify_ref(c)
    a = c%4
    b = a+c÷4
    if b>3
        b = b%3
    end
    return [a,b]
end

"""

    h_conformity!(mesh::HPMesh)

Performs the marking of previously un-marked triangles (and corresponding edges) in order to avoid hanging nodes. This function must be run after the marking of the unsuitable triangles, by `mark_triangles`.
"""
function h_conformity!(mesh::HPMesh)
    (;ET,rE,rT,leT) = mesh 
    nm0  = 0
    nm   = sum(rE)
    while nm != nm0
        nm0 = nm 
        for k in eachindex(rT)
            sum(rE[ET[:,k]])>0 && (rE[ET[leT[k],k]]=true)
        end
        nm = sum(rE)
    end
    for k in eachindex(rT)
        marked = findall(rE[ET[:,k]])
        if !isempty(marked)
            s = length(marked)
            if s == 1
                rT[k] = leT[k]
            elseif  s == 2
                rT[k] = codify_ref(leT[k],setdiff(marked,leT[k])[1])
            else
                rT[k] = 4
            end
        end
    end
end



"""

    find_edge(e,E) 

returns the index of the column of `E` formed by `e`, or `size(E,2)+1` if `e` is not in any column of `E`.
"""
function find_edge(e,E) 
    i = findfirst(z->issubset(e,z),eachcol(E))
    i === nothing ? size(E,2)+1 : i
end


"""

    green_ref!(mesh::HPMesh,i::RefInd,omesh::HPMesh,k)

performs a green refinement of triangle `k` in mesh `omesh`, saving the output in `mesh`. `i` keeps track of the indices of the refinement process. 
"""
function green_ref!(mesh::HPMesh,i::RefInd,omesh::HPMesh,k)
    println("Entering Green, triangle $k")
    (;P,T,E,ET,pE,leT) = mesh
    e             = omesh.ET[omesh.rT[k],k] #the edge to be refined
    oe            = setdiff(omesh.ET[:,k],e) #the other edges, that should be copied
    i.nodes[1:2] .= omesh.E[:,e] #indices of nodes in the refined edge
    i.nodes[3]    = setdiff(omesh.T[:,k],omesh.E[:,e])[1] # index of node opposite to the refined edge
    # We need to know which of the other edges correspond to each of the nodes of the refined edge. 
    #  the following line makes sure that oe[1] ∋ i.nodes[1].
    i.nodes[2] ∈ omesh.E[:,oe[1]] ? reverse!(oe) : nothing  
    # Now, it is possible that e has already been visited and refined. In that case, i.EE[e]!=0 
    #  and we need to recover the information about the new edges.
    #  For this, we assume that the new edges where created following the order of e. In other words, 
    #  omesh.E[:,e] = [i.nodes[1],i.nodes[2]], so the first new edge contains i.nodes[1] and the new node i.nodes[4],
    #  and the second new edge contains i.nodes[2] and i.nodes[4].
    if i.EE[e] > 0 
            println("edge has been previously seen...")
        i.nodes[4] = i.PE[e]
        i.edges[1] = i.EE[e]
        i.edges[2] = i.EE[e]+1
        i.edges[3] = i.PTE[3]
        i.PTE[3]  += 1
        # We add the new edge (the interior one):
        E[:,i.edges[3]] = [i.nodes[2],i.nodes[4]]
    # If e has not been refined yet, we refine it:
    else
        P[:,i.PTE[1]]     = sum(mesh.P[:,omesh.E[:,e]],dims=2)/2 #new node
        i.nodes[4]   = i.PTE[1]
        i.PE[e]      = i.PTE[1]
        i.PTE[1]    += 1
        i.EE[e]      = i.PTE[3]
        i.edges[1:3].= [i.PTE[3],i.PTE[3]+1,i.PTE[3]+2] 
        i.PTE[3]    += 3
        #We add the new edges, where the edge e is divided following the original order: e[1]->new->e[2].
        E[:,i.edges[1:3]] = i.nodes[[1 4 2;
                                     4 2 4]]
    end
    # Now, we check if the other edges have already been visited. If they were, we recover the corresponding index
    #  from i.EE, taking into accout that the index can be stored as negative if the edge was save in reverse order. 
    for (j,edge) in enumerate(oe)
        if i.EE[edge]!=0
            i.edges[3+j] = i.EE[edge]
        else 
            i.edges[3+j]   = i.PTE[3]
            E[:,i.PTE[3]]  = omesh.E[:,edge]
            i.EE[edge]     = i.PTE[3]
            i.PTE[3]      += 1
        end
    end
            
    T[:,i.PTE[2]:i.PTE[2]+1]  .= i.nodes[[3 2; 
                                          1 3;
                                          4 4]]            
    ET[:,i.PTE[2]:i.PTE[2]+1] .= i.edges[[5 4;
                                          1 3;
                                          3 2]]
    leT[i.PTE[2]:i.PTE[2]+1]  .= ones(2)
    i.PTE[2]  += 2
    pE[i.edges[1:2]] .= omesh.pE[e]
    pE[i.edges[4]] = omesh.pE[oe[1]]
    pE[i.edges[5]] = omesh.pE[oe[2]]
    pE[i.edges[3]] = max(1,abs(pE[i.edges[5]]-pE[i.edges[1]]),abs(pE[i.edges[4]]-pE[i.edges[2]]))
end
    

"""

    blue_ref!(mesh::HPMesh,i::RefInd,omesh::HPMesh,k)

performs a blue refinement of triangle `k` in mesh `omesh`, saving the output in `mesh`. `i` keeps trak of the indices of the refinement process. 
"""
function blue_ref!(mesh::HPMesh,i::RefInd,omesh::HPMesh,k)
    println("Entering Blue, triangle $k")
    (;P,T,E,ET,pE,leT) = mesh
    et      = omesh.ET[decodify_ref(omesh.rT[k]),k]
    edges   = omesh.E[:,et]
    # nodes are stored in order: first node of main edge, shared node between edges to be refined, third node. 
    i.nodes[1] = setdiff(edges[:,1],edges[:,2])[1]
    i.nodes[2] = setdiff(edges[:,1],i.nodes[1])[1]
    i.nodes[3] = setdiff(edges[:,2],i.nodes[1:2])[1]
    oe = setdiff(omesh.ET[:,k],et)[1]

    for (j,e) in enumerate(et)
        if i.EE[e] > 0 
            i.nodes[3+j] = i.PE[e]
            if i.nodes[j:j+1] == edges[:,j]
                i.edges[2j-1:2j] .= i.EE[e]:i.EE[e]+1
            else 
                i.edges[2j-1:2j] .= [i.EE[e]+1,i.EE[e]]    
            end
        else
            P[:,i.PTE[1]]    = sum(mesh.P[:,omesh.E[:,e]],dims=2)/2 #new node
            i.PE[e]       = i.PTE[1]
            i.nodes[3+j]  = i.PTE[1]
            i.PTE[1]     += 1
            i.EE[e]       = i.PTE[3]
            if i.nodes[j:j+1] == omesh.E[:,e]
                i.edges[2j-1:2j]     .= i.PTE[3]:i.PTE[3]+1
                E[:,i.edges[2j-1:2j]] = i.nodes[[j 3+j;
                                                 3+j j+1]]
            else 
                i.edges[2j-1:2j]     .= i.PTE[3]+1:-1:i.PTE[3]
                E[:,i.edges[2j-1:2j]] = i.nodes[[j+1 3+j;
                                                 3+j  j]]
            end
            i.PTE[3]  += 2
        end
    end
    if i.EE[oe] !=0 
        i.edges[5] = i.EE[oe]
    else 
        i.edges[5] = i.PTE[3]
        i.EE[oe]   = i.PTE[3]  
        i.PTE[3]  += 1
        E[:,i.edges[5]] = omesh.E[:,oe]
    end
    i.edges[6:7]     .= i.PTE[3]:i.PTE[3]+1 
    i.PTE[3]         += 2
    E[:,i.edges[6:7]] = i.nodes[[3 4;
                                 4 5]]    
    T[:,i.PTE[2]:i.PTE[2]+2]  = i.nodes[[3 4 3; 
                                         1 2 4;
                                         4 5 5]]            
    ET[:,i.PTE[2]:i.PTE[2]+2] = i.edges[[5 2 6;
                                         1 3 7;
                                         6 7 4]]
    leT[i.PTE[2]:i.PTE[2]+2]   .= ones(3)
    pE[i.edges[1:2]] .= omesh.pE[omesh.ET[et[1]]]
    pE[i.edges[3:4]] .= omesh.pE[omesh.ET[et[2]]]
    pE[i.edges[5]]    = omesh.pE[oe]
    pE[i.edges[6]]    = max(1,abs(pE[i.edges[1]]-pE[i.edges[5]]),abs(pE[i.edges[2]]-pE[i.edges[4]]))
    pE[i.edges[7]]    = max(1,abs(pE[i.edges[2]]-pE[i.edges[3]]),abs(pE[i.edges[4]]-pE[i.edges[6]]))
    i.PTE[2]  += 3

end


"""

    red_ref!(mesh::HPMesh,i::RefInd,omesh::HPMesh,k)

performs a red refinement of triangle `k` in mesh `omesh`, saving the output in `mesh`. `i` keeps track of the indices of the refinement process, through `mark_triangles` and `h_conformity`. 
"""
function red_ref!(mesh::HPMesh,i::RefInd,omesh::HPMesh,k)
    println("Entering Red, triangle $k")
    (;P,T,E,ET,pE,leT) = mesh
    # et are the indices of the edges beginning by the longest one. 
    et           = omesh.ET[mod1.(omesh.leT[k]:omesh.leT[k]+2,3),k]
    i.nodes[1]   = setdiff(omesh.E[:,et[1]],omesh.E[:,et[2]])[1]
    i.nodes[2]   = setdiff(omesh.E[:,et[1]],i.nodes[1])[1]
    i.nodes[3]   = setdiff(omesh.E[:,et[2]],i.nodes[1:2])[1]
    for (j,e) in enumerate(et)
        e = et[j]
        if i.EE[e] > 0 
            i.nodes[3+j] = i.PE[e]
            if [i.nodes[j],i.nodes[mod1(j+1,3)]] == omesh.E[:,e]
                i.edges[2j-1:2j] .= i.EE[e]:i.EE[e]+1
            else 
                i.edges[2j-1:2j] .= [i.EE[e]+1,i.EE[e]]    
            end
        else
            P[:,i.PTE[1]]= sum(mesh.P[:,omesh.E[:,e]],dims=2)/2 #new node
            i.PE[e]      = i.PTE[1]
            i.nodes[3+j] = i.PTE[1]
            i.PTE[1]    += 1
            i.EE[e]      = i.PTE[3]
            if [i.nodes[j],i.nodes[mod1(j+1,3)]] == omesh.E[:,e]
                i.edges[2j-1:2j] .= i.PTE[3]:i.PTE[3]+1
                E[:,i.edges[2j-1:2j]] = i.nodes[[ j  3+j;
                                                 3+j mod1(j+1,3)]]
            else 
                i.edges[2j-1:2j] .= [i.PTE[3]+1,i.PTE[3]]
                E[:,i.edges[2j-1:2j]] = i.nodes[[mod1(j+1,3) 3+j; 
                                                 3+j j]]
            end
            i.PTE[3]  += 2
        end
    end
    i.edges[7:9] .= i.PTE[3]:i.PTE[3]+2
    i.PTE[3]     += 3
    E[:,i.edges[7:9]] = i.nodes[[6 4 5;
                                 4 5 6]]
    T[:,i.PTE[2]:i.PTE[2]+3]  = i.nodes[[1 4 6 5; 
                                         4 2 5 6;
                                         6 5 3 4]]            
    ET[:,i.PTE[2]:i.PTE[2]+3] = i.edges[[1 2 9 9;
                                         7 3 4 7;
                                         6 8 5 8]]
    
    leT[i.PTE[2]:i.PTE[2]+3]   .= ones(4)
    pE[i.edges[1:2]] .= omesh.pE[et[1]]
    pE[i.edges[3:4]] .= omesh.pE[et[2]]
    pE[i.edges[5:6]] .= omesh.pE[et[3]]
    pE[i.edges[7]]    = max(1,abs(pE[i.edges[1]]-pE[i.edges[6]]))
    pE[i.edges[8]]    = max(1,abs(pE[i.edges[2]]-pE[i.edges[3]]))
    pE[i.edges[9]]    = max(1,abs(pE[i.edges[4]]-pE[i.edges[5]]))
    i.PTE[2]  += 4
    println("exit red, i.EE:")
    display(Int.(i.EE))
end

function no_ref!(mesh::HPMesh,i::RefInd,omesh::HPMesh,k)
    (;T,E,ET,leT,pE) = mesh
    et           = omesh.ET[:,k]
    i.nodes[1:3] = omesh.T[:,k]
    for j in 1:3
        e = et[j]
        if i.EE[e] != 0 
            i.edges[j] = i.EE[e]
        else
            i.EE[e]    = i.PTE[3]
            i.edges[j] = i.PTE[3]
            i.PTE[3]      += 1
            pE[i.edges[j]] = omesh.pE[e]
        end
    end
    E[:,i.edges[1:3]] .= omesh.E[:,et] 
    T[:,i.PTE[2]]          .= omesh.T[:,k]
    ET[:,i.PTE[2]]         .= i.edges[1:3]
    leT[i.PTE[2]]           = omesh.leT[k]
    i.PTE[2]  += 1
end
"""

    refine_h(omesh::HPMesh)

performs the `h` refinement of `omesh`, returning a new mesh. `omesh.rE` and `omesh.rT` should be properly loaded with the refinement information  
"""
function refine_h(omesh::HPMesh)
    red          = findall(omesh.rT.==4)
    blue         = findall(omesh.rT.>4)
    green        = findall(1 .<=omesh.rT.<=3)    
    
    n_new_tris   = size(omesh.T,2) + 3length(red) + 2length(blue) + length(green)
    n_new_points = size(omesh.P,2) + sum(omesh.rE)
    n_new_edges  = length(omesh.rE) + sum(omesh.rE) + 3length(red) + 2length(blue) + length(green)
    mesh   = HPMesh(zeros(2,n_new_points), #P
                    zeros(UInt32,3,n_new_tris),#T,
                    zeros(UInt32,2,n_new_edges),#E,
                    zeros(UInt32,3,n_new_tris),#ET,
                    zeros(UInt8,n_new_points) ,#Pm,
                    zeros(UInt8,n_new_edges),#Em,
                    zeros(UInt8,n_new_edges) ,#pE,
                    zeros(UInt8,3,n_new_tris),#pT,
                    zeros(UInt8,n_new_tris) ,#rT,
                    zeros(Bool,n_new_edges) ,#rE,
                    ones(UInt8,n_new_tris) , #leT 
                    nothing,
                    nothing)
    mesh.P[:,1:size(omesh.P,2)] = omesh.P
    println(typeof(mesh)==typeof(omesh))
    i   = RefInd(MVector{3}{UInt32}(UInt32[size(omesh.P,2)+1,1,1]),
                 MVector{length(omesh.rE)}{UInt32}(zeros(UInt32,length(omesh.rE))),
                 MVector{length(omesh.rE)}{UInt32}(zeros(UInt32,length(omesh.rE))),
                 MVector{6}{UInt32}(zeros(UInt32,6)), 
                 MVector{9}{UInt32}(zeros(UInt32,9)), )
    for k in eachindex(omesh.rT)
        if omesh.rT[k] == 0

            no_ref!(mesh,i,omesh,k)
        elseif 1<=omesh.rT[k]<=3
            green_ref!(mesh,i,omesh,k)
        elseif omesh.rT[k]>4
            blue_ref!(mesh,i,omesh,k)
        elseif omesh.rT[k] == 4
            red_ref!(mesh,i,omesh,k)
        end
    end
    return mesh    
end

function refine_p!(mesh::HPMesh)
    #Esta función debe recorrer los elementos marcados para refinamiento p que no fueron previamente refinados h. Para ello, el refinado h debe recibir el marcado p y pasar a false los elementos que se parten. O bien este proceso debe realizarse previamente por fuera del refinado, pero luego de h_conformity. 
    false
end


function check_conformity(pt::Vector{T}) where T<:Unsigned
    sum(pt) - maximum(pt) ≥ maximum(pt) 
end


function p_conformity!(mesh::HPMesh,k::U,d) where U<:Unsigned    
    (;P,T,E,pE) = mesh
    pt  = view(mesh.pT,:,k)
    et  = view(mesh.ET,:,k)
    i,j = partialsortperm(pt,1:2)
    out = false
    if check_comformity(pT[:,k])
        out = true
    else
        if d>0
            pt[i]     += 1
            pE[et[i]] += 1
            k₁,k₂      = 
            if p_conformity!(mesh,k₁,d-1)
                out = true
            else
                pt[i]     -= 1
                pE[et[i]] -= 1
                pt[j]     += 1
                pE[et[j]] += 1
                if p_conformity!(mesh,k₂,d-1)
                    out = true
                else
                    pt[j] -= 1
                    pE[et[j]] -= 1
                end
            end
        end
    end
    out
end


function jump()
    
end

