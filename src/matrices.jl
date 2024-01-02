""" 
    `MatrixDict` and `ArrayDict` are unused, yet. The idea is to store the local mass and stiff matrices for each combination of values of `p`.
""" 

const MatrixDict = Dict{NTuple{3,UInt8},SMatrix{Float64}}
const ArrayDict  = Dict{NTuple{3,UInt8},SArray}



function degrees_of_freedom_by_edge(mesh::MeshHP{F,I,P,Bool}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;edges) = mesh 
    by_edge  = similar(mesh.edges,Vector{I})
    i        = size(mesh.points,2)+1
    for edge in edges
        nod  = nodes(edge)
        med  = collect(i:i+degree(edge)-2)
        set!(by_edge,edge,[nod[1],med...,nod[2]])
        i   += degree(edge)-1
    end
    return by_edge,i
end

function degrees_of_freedom(mesh::MeshHP{F,I,P,Bool}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;edges,triangles) = mesh
    by_edge,i = degrees_of_freedom_by_edge(mesh)
    dof     = similar(mesh.triangles,Vector{I})
    for t in triangles
        edges_t      = get_edges(t)
        true_edges_t = [edges[e] for e in edges_t]
        pind         = psortperm(degree(e) for e in true_edges_t)        
        dof[t]       = zeros(I,compute_dimension((degree(e) for e in true_edges_t)...))
        j = 1
        @inbounds for i in pind
            newdof = by_edge[edges_t[i]]
            if edges_t[i] == nodes(true_edges_t[i])
                dof[t][j:j+length(newdof)-1] .= newdof[1:end-1]
            else
                dof[t][j:j+length(newdof)-1] .= reverse(newdof[2:end])
            end 
            j += length(newdof)
        end
        dof[t][j:end] = j:j+length(dof[t])
    end
    return dof
end

function transform_matrix!(A,vert) 
    @views A[:,1] .= vert[:,2]-vert[:,1] 
    @views A[:,2] .= vert[:,3]-vert[:,1]
end


"""
    compute_hat_mass(B::Basis)

computes the mass matrix in the reference triangle fot the basis `B`. 
"""
function compute_hat_mass(B::Basis)
    (;nodes,b,C) = B 
    sch = grundmann_moeller(Float64,Val(2),2length(nodes)+1)
    n   = length(b)
    M   = zeros(n,n)
    T   = [[-1,-1],[1,-1],[-1,1]]
    @inbounds for j in 1:n, i in 1:j
        M[i,j] = integrate(x->b[i](x)*b[j](x),sch,T)
    end 
    return SMatrix{n,n}(C'*Symmetric(M)*C)
end;



"""
    compute_hat_stiff(B::Basis)

computes the stiff matrix in the reference triangle fot the basis `B`. 
"""
function compute_hat_stiff(base::Basis)
    (;nodes,∇b,C) = base
    sch = grundmann_moeller(Float64,Val(2),2*length(nodes)+1)
    n   = length(∇b)
    S   = zeros(n,n,4)
    T   = [[-1,-1],[1,-1],[-1,1]]
    @inbounds for j in 1:n, i in 1:n
        S[i,j,1] = integrate(x->∇b[i](x)[1]*∇b[j](x)[1],sch,T)
        S[i,j,2] = integrate(x->∇b[i](x)[2]*∇b[j](x)[1],sch,T)
        S[i,j,3] = integrate(x->∇b[i](x)[1]*∇b[j](x)[2],sch,T)
        S[i,j,4] = integrate(x->∇b[i](x)[2]*∇b[j](x)[2],sch,T)
    end
    #@inbounds for k in 1:4
    #    S[:,:,k] = Symmetric(S[:,:,k])
    #end
    SArray{Tuple{n,n,4},Float64}(S)
end;



function is_in_0(v,V;tol=1e-14)
    i     = 1
    while i≤size(V,2)
        norm(v-V[:,i])<tol && return i
        i +=1 
    end
    return 0
end


function compute_mass(mesh::MeshHP,dof,MD::MatrixDict)
    (;points,triangles,edges) = mesh
    ℓ = sum(length.(dofs))
    I = Vector{Int64}(undef,ℓ)
    J = Vector{Int64}(undef,ℓ)
    V = Vector{Float64}(undef,ℓ)
    Aₜ= MMatrix{2,2}(zeros(2,2))
    r = 1
    @inbounds for t in triangles
        dofT = dof[t] 
        ℓₜ   = length(dofT)
        trasform_matrix!(Aₜ,view(points,:,nodes(t)))
        dAₜ  = abs(0.25det(Aₜ))
        p    = []# HAY QUE CALCULAR N PARA TRES p DISTINTOS Y VER SI NECESITO LOS DOFS COMPLETOS O Qué ..
        n    = Int((2p[2]-p[1]+2)*(p[1]+1)/2)
        v    = dAₜ*MD[p] 
        i = repeat(dofT,n)
        j = repeat(dofT,inner=n)
        I[r:r+ℓₜ^2-1] = i
        J[r:r+ℓₜ^2-1] = j
        V[r:r+ℓₜ^2-1] = v
        r += ℓₜ^2
    end
    sparse(I,J,V)
end



"""
    compute_stiff_std(mesh::MeshHP,BD::BasisDict)

computes the stiffness matrix corresponing to `mesh` following a standard approach i.e,: integrating on each triangle and not using the hat matrix.  
"""

function compute_stiff_std(mesh::MeshHP,BD::BasisDict)
    (;P,T,pT,pE,D,Da) = mesh
    ℓ = sum(x->length(x)^2,Da)
    I = Vector{Int64}(undef,ℓ)
    J = Vector{Int64}(undef,ℓ)
    V = Vector{Float64}(undef,ℓ)
    Aₜ= MMatrix{2,2}(zeros(2,2))
    iAₜ = MMatrix{2,2}(zeros(2,2))
    r = 1
    pmax = maximum(pE)
    sch = grundmann_moeller(Float64,Val(2),4pmax+1)
    t̂ = [[-1.,-1.],[1.,-1.],[-1.,1.]]
    @inbounds for k in eachindex(eachcol(T))
        p    = NTuple{3,Int8}(pT[:,k])
        bd   = BD[p]
        t    = view(T,:,k)
        dofT = Da[k]
        ℓₜ   = length(dofT)
        Aₜ  .= transform_matrix(P,t)
        iAₜ .= inv(Aₜ)
        dAₜ  = abs(det(Aₜ))
        C    = bd.C
        ∇b   = bd.∇b
        v    = zeros(ℓₜ,ℓₜ)
        for j in 1:ℓₜ, i in 1:j
            v[i,j] = dAₜ*integrate(x->(iAₜ'*∇b[i](x))⋅(iAₜ'*∇b[j](x)),sch,t̂)
        end
        v   = C'*Symmetric(v)*C 
        i = repeat(dofT,ℓₜ)
        j = repeat(dofT,inner=ℓₜ)
        I[r:r+ℓₜ^2-1] = i
        J[r:r+ℓₜ^2-1] = j
        V[r:r+ℓₜ^2-1] = v[:]
        r += ℓₜ^2
    end
    sparse(I,J,V)
end


""" 
    compute_rhs(mesh,BD,f)

computes the right-hand side, ∫ϕf corresponing to `mesh` if data `f`. 
"""
function compute_rhs(mesh,BD,f)
    (;P,T,D,Da,pE,pT) = mesh
    ℓ = sum(length,Da)
    pmax = maximum(pE)
    #Constructors
    I = Vector{Int64}(undef,ℓ)
    V = Vector{Float64}(undef,ℓ)
    Aₜ   = MMatrix{2,2}(zeros(2,2))
    bₜ   = MVector{2}(zeros(2))
    t̂    = [[-1. -1.],[1. -1.],[-1. 1.]]
    sch  = grundmann_moeller(Float64,Val(2),4pmax+1)
    r    = 1
    for k in eachindex(eachcol(T))
        p    = NTuple{3,Int8}(pT[:,k])
        t    = view(T,:,k)
        Aₜ  .= transform_matrix(P,t)
        dAₜ  = abs(det(Aₜ))
        bₜ  .= transform_term(P,t)
        b    = view(BD[p].b,:)
        C    = view(BD[p].C,:,:)
        ℓₜ   = length(b)
        F    = zeros(ℓₜ)
        for i in eachindex(F)
            F[i] = dAₜ*integrate(x->b[i](x)*f(Aₜ*x+bₜ),sch,t̂)
        end
        v    = C'*F
        I[r:r+ℓₜ-1] .= Da[k]
        V[r:r+ℓₜ-1] .= v
        r += ℓₜ
    end
    sparsevec(I,V,length(D))
end


""" 

    hat_local_mass(p₁[,p₂])

Returns the local mass matrix on the reference triangle, with vertices [-1,-1], [1,-1] and [-1,1].

""" 
function hat_local_mass(p₁,p₂)
    sch = grundmann_moeller(Float64,Val(2),2p₂+1)
    B   = standard_basis(p₁,p₂)
    n   = length(B)
    M   = zeros(n,n)
    T̂   = [[-1,-1],[1,-1],[-1,1]]
    for i in 1:n
        for j in i:n
            M[i,j] = integrate(x->B[i](x)*B[j](x),sch,T̂)
        end
    end
    nodes = boundary_nodes(p₁,p₂)
    F     = matrix_F(B,nodes)
    C     = matrix_C(F)
    return C'*Symmetric(M)*C
end;
hat_local_mass(p₁) = hat_local_mass(p₁,p₁);
