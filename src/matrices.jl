""" 
    `MatrixDict` and `ArrayDict` are unused, yet. The idea is to store the local mass and stiff matrices for each combination of values of `p`.
using Base: io_has_tvar_name
""" 

const MatrixDict{I} = Dictionary{DegTuple{I},SMatrix{Float64}}
const ArrayDict{I}  = Dictionary{DegTuple{I},SArray}



function degrees_of_freedom_by_edge(mesh::MeshHP{F,I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;edgelist) = mesh 
    by_edge  = similar(mesh.edgelist,Vector{I})
     i        = size(mesh.points,2)+1
    for edge in edges(edgelist)
        nod  = nodes(edge)
        med  = collect(i:i+degree(edgelist[edge])-2)
        set!(by_edge,edge,[nod[1],med...,nod[2]])
         i   += degree(edgelist[edge])-1
    end
    return by_edge
end

function degrees_of_freedom(mesh::MeshHP{F,I,P}) where {F,I,P}
    by_edge = degrees_of_freedom_by_edge(mesh)
    degrees_of_freedom(mesh,by_edge)
end

function degrees_of_freedom(mesh::MeshHP{F,I,P},by_edge::Dictionary{EdgeHP{I},Vector{I}}) where {F,I,P}
    (;edgelist,trilist) = mesh
    dof     = similar(mesh.trilist,Vector{I})
    k       = maximum(maximum.(by_edge))+1 #first non-edge dof
    for t in triangles(trilist)
        p,t_edges = pedges(t,mesh)
        dof[t]  = zeros(I,compute_dimension(p))
        j = 1 #counter of dof in current triangle
        @inbounds for i in 1:3
            newdof = by_edge[t_edges[i]]
            if  same_order(t_edges[i],edgelist)
                dof[t][j:j+length(newdof)-2] .= newdof[1:end-1]
            else
                dof[t][j:j+length(newdof)-2] .= reverse(newdof[2:end])
            end 
            j += length(newdof)-1
        end
        dof[t][j:end] = k:k+(length(dof[t])-j)
        k += length(dof[t])-j + 1
    end
    return dof
end

function boundary_dof(mesh::MeshHP{F,I,P},by_edge) where {F,I,P}
    (;edgelist) = mesh
    v = I[]
    for e in edges(edgelist)
        if marker(e)==1
            push!(v,by_edge[e]...)
        end
    end
    return v
end

function transform_matrix!(A,vert) 
    @views A[:,1] .= 0.5(vert[:,1]-vert[:,2]) 
    @views A[:,2] .= 0.5(vert[:,3]-vert[:,2])
end

function transform_term!(b,vert)
    @views b .= 0.5(vert[:,1]+vert[:,3])
end

"""
    compute_hat_mass(B::Basis)

computes the mass matrix in the reference triangle fot the basis `B`. 
"""
function compute_hat_mass(B::Basis)
    (;dim,b,C) = B 
    sch = grundmann_moeller(Float64,Val(2),2dim+1)
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
    (;dim,∇b) = base
    sch = grundmann_moeller(Float64,Val(2),4*dim+1)
    S   = zeros(dim,dim,4)
    t̂   = [[-1,1],[-1,-1],[1,-1]]
    @inbounds for j in 1:dim, i in 1:j
        S[i,j,1] = integrate(x->∇b[i](x)[1]*∇b[j](x)[1],sch,t̂)
        S[i,j,2] = integrate(x->∇b[i](x)[2]*∇b[j](x)[1],sch,t̂)
        S[i,j,3] = integrate(x->∇b[i](x)[1]*∇b[j](x)[2],sch,t̂)
        S[i,j,4] = integrate(x->∇b[i](x)[2]*∇b[j](x)[2],sch,t̂)
    end
    @inbounds for k in 1:4
        S[:,:,k] = Symmetric(S[:,:,k])
    end
    SArray{Tuple{dim,dim,4},Float64}(S)
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
    (;points,trilist) = mesh
    ℓ = maximum(maximum.(dof))
    I = Vector{Int64}(undef,ℓ)
    J = Vector{Int64}(undef,ℓ)
    V = Vector{Float64}(undef,ℓ)
    Aₜ= MMatrix{2,2}(zeros(2,2))
    r = 1
    @inbounds for t in trilist
        dofT = dof[t] 
        ℓₜ   = length(dofT)
        p,pnod  = pnodes(t,mesh)
        transform_matrix!(Aₜ,view(points,:,pnod))
        dAₜ  = abs(det(Aₜ))
        n    = compute_dimension(p)
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



function stiff(mesh::MeshHP,dof,AD::ArrayDict,bd::BasisDict)
    (;points,trilist) = mesh
    ℓ = sum(x->length(x)^2,dof)
    #ℓ = maximum(maximum.(dof))
    I = Vector{Int64}(undef,ℓ)
    J = Vector{Int64}(undef,ℓ)
    V = Vector{Float64}(undef,ℓ)
    Aₜ= MMatrix{2,2}(zeros(2,2))
    iAₜ= MMatrix{2,2}(zeros(2,2))
    # iAₜ² = MMatrix{2,2}(zeros(2,2))
    r = 1
    @inbounds for t in triangles(trilist)
        dofT = dof[t] 
        p,pnod  = pnodes(t,mesh)
        (;dim,C) = bd[p]
        transform_matrix!(Aₜ,view(points,:,pnod))
        iAₜ .= inv(Aₜ)
        iAₜ .= iAₜ'*iAₜ
        z = vec(iAₜ)
        dAₜ  = abs(det(Aₜ))
        v    = zeros(dim,dim)
        S = AD[p]
        for j in 1:dim, i in 1:j
            v[i,j] = S[i,j,:]⋅z
        end
        v = dAₜ*C'*Symmetric(v)*C
        i = repeat(dofT,dim)
        j = repeat(dofT,inner=dim)
        I[r:r+dim^2-1] = i
        J[r:r+dim^2-1] = j
        V[r:r+dim^2-1] = v
        r += dim^2
    end
    sparse(I,J,V)
end
"""
    stiff_std(mesh::MeshHP,BD::BasisDict)

computes the stiffness matrix corresponing to `mesh` following a standard approach i.e,: integrating on each triangle and not using the hat matrix.  
"""

function stiff_std(mesh::MeshHP,dof,BD::BasisDict)
    (;points,trilist,edgelist) = mesh
    ℓ = sum(x->length(x)^2,dof)
    I = Vector{Int64}(undef,ℓ)
    J = Vector{Int64}(undef,ℓ)
    V = Vector{Float64}(undef,ℓ)
    Aₜ= MMatrix{2,2}(zeros(2,2))
    iAₜ = MMatrix{2,2}(zeros(2,2))
    r = 1
    pmax = maximum(degree.(edgelist))
    sch = grundmann_moeller(Float64,Val(2),4pmax+1)
    t̂ = [[-1.,1.],[-1.,-1.],[1.,-1.]]
    @inbounds for t in triangles(trilist)
        p,pnod  = pnodes(t,mesh)
        (;dim,∇b,C)  = BD[p]
        dofT = dof[t]
        transform_matrix!(Aₜ,view(points,:,pnod))
        iAₜ .= inv(Aₜ)
        dAₜ  = abs(det(Aₜ))
        v    = zeros(dim,dim)
        for j in 1:dim, i in 1:j
            v[i,j] = dAₜ*integrate(x->(iAₜ*∇b[i](x))⋅(iAₜ*∇b[j](x)),sch,t̂)
        end
        v   = C'*Symmetric(v)*C 
        I[r:r+dim^2-1] = repeat(dofT,dim)
        J[r:r+dim^2-1] = repeat(dofT,inner=dim)
        V[r:r+dim^2-1] = v[:]
        r += dim^2
    end
    sparse(I,J,V)
end

# AL HACER transform_matrix! HAY QUE CONTAR LOS NODES ORDENADOS SEGUN p. 

""" 
    compute_rhs(mesh,BD,f)

computes the right-hand side, ∫ϕf corresponing to `mesh` if data `f`. 
"""
function rhs(mesh,dof,BD,f)
    (;points,edgelist) = mesh
    ℓ    = sum(length,dof)
    pmax = maximum(degree.(edgelist))
    #Constructors
    I   = Vector{Int64}(undef,ℓ)
    V   = Vector{Float64}(undef,ℓ)
    Aₜ  = MMatrix{2,2}(zeros(2,2))
    bₜ  = MVector{2}(zeros(2))
    t̂ = [[-1.,1.],[-1.,-1.],[1.,-1.]]
    sch = grundmann_moeller(Float64,Val(2),4pmax+1)
    r   = 1
    @inbounds for t in triangles(mesh)
        p,pnod  = pnodes(t,mesh)
        transform_matrix!(Aₜ,view(points,:,pnod))
        transform_term!(bₜ,view(points,:,pnod))
        dAₜ  = abs(det(Aₜ))
        (;b,C) = BD[p] 
        ℓₜ   = length(b)
        F    = zeros(ℓₜ)
        for i in eachindex(F)
            F[i] = dAₜ*integrate(x->b[i](x)*f(Aₜ*x+bₜ),sch,t̂)
        end
        v    = C'*F
        I[r:r+ℓₜ-1] .= dof[t]
        V[r:r+ℓₜ-1] .= v
        r += ℓₜ
    end
    Vector(sparsevec(I,V,maximum(maximum.(dof))))
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
    t̂ = [[-1.,1.],[-1.,-1.],[1.,-1.]]
    for i in 1:n
        for j in i:n
            M[i,j] = integrate(x->B[i](x)*B[j](x),sch,t̂)
        end
    end
    nodes = boundary_nodes(p₁,p₂)
    F     = matrix_F(B,nodes)
    C     = matrix_C(F)
    return C'*Symmetric(M)*C
end;
hat_local_mass(p₁) = hat_local_mass(p₁,p₁);

function basisdict(mesh::MeshHP{F,I,P}) where {F,I,P}
    bd = BasisDict{P}()
    for t in triangles(mesh)
        p,_ = pnodes(t,mesh)
        isin,_ = gettoken(bd,p)
        if !isin
            set!(bd,p,Basis(p))
        end
    end
    return bd
end
function stiffdict(mesh::MeshHP{F,I,P},bd::BasisDict) where {F,I,P}
    ad = ArrayDict{P}()
    for t in triangles(mesh)
        p,_ = pnodes(t,mesh)
        isin,_ = gettoken(ad,p)
        if !isin
            set!(ad,p,compute_hat_stiff(bd[p]))
        end
    end
    return ad
end

function massdict(mesh::MeshHP{F,I,P},bd::BasisDict) where {F,I,P}
    md = MatrixDict{P}()
    for t in triangles(mesh)
        p,_ = pnodes(t,mesh)
        isin,_ = gettoken(md,p)
        if !isin
            set!(md,p,compute_hat_mass(bd[p]))
        end
    end
    return md
end
