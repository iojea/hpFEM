""" 
    `MatrixDict` and `ArrayDict` are unused, yet. The idea is to store the local mass and stiff matrices for each combination of values of `p`.
"""

const MatrixDict{P} = Dictionary{DegTuple{P},SMatrix{Float64}}
const BiMatrixDict{P} = Dictionary{DegTuple{P},NTuple{2,SMatrix{Float64}}}
const ArrayDict{P} = Dictionary{DegTuple{P},SArray}
const DofDict{P} = Dictionary{TriangleHP{P},Vector{P}}
const EdgeDofDict{P} = Dictionary{EdgeHP{P},Vector{P}}

"""
    degrees_of_freedom_by_edge(mesh::MeshHP{F,I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}

Creates a dictionary (from `Dictionaries.jl`) where the keys are the edges of `mesh` and the values are vectors with indices corresponding to the nodal degrees of freedom. 
"""
function degrees_of_freedom_by_edge(
    mesh::MeshHP{F,I,P},
) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (; edgelist) = mesh
    by_edge = similar(mesh.edgelist, Vector{I})
    i = size(mesh.points, 2) + 1
    for edge in edges(edgelist)
        med = collect(i:i+degree(edgelist[edge])-2)
        set!(by_edge, edge, [edge[1], med..., edge[2]])
        i += degree(edgelist[edge]) - 1
    end
    return by_edge
end

"""
    degrees_of_freedom(mesh::MeshHP{F,I,P}) where {F,I,P}

Creates a dictionary (from `Dictionaries.jl`) where the keys are the triangles of the mesh and the values are vectores storing the indices of the corresponding degrees of freedom. 

Internally, `degrees_of_freedom_by_edge` is called in order to obtain the nodal degrees of freedom, and then the degrees of freedom corresponding to bubble functions are computed. 

If a dictionary of degrees of freedom by edge has already been computed, it is recommended to run: 

    degrees_of_freedom(mesh::MeshHP{F,I,P},by_edge::Dictionary{EdgeHP{I},Vector{I}}) where {F,I,P}
"""
function degrees_of_freedom(mesh::MeshHP{F,I,P}) where {F,I,P}
    by_edge = degrees_of_freedom_by_edge(mesh)
    degrees_of_freedom(mesh, by_edge)
end
function degrees_of_freedom(
    mesh::MeshHP{F,I,P},
    by_edge::Dictionary{EdgeHP{I},Vector{I}},
) where {F,I,P}
    (; edgelist, trilist) = mesh
    dof = similar(trilist, Vector{I})
    k = maximum(maximum.(by_edge)) + 1 #first non-edge dof
    for t in triangles(trilist)
        p, t_edges = pedges(t, mesh)
        dof[t] = zeros(I, compute_dimension(p))
        j = 1 #counter of dof in current triangle
        @inbounds for i = 1:3
            newdof = by_edge[t_edges[i]]
            if same_order(t_edges[i], edgelist)
                dof[t][j:j+length(newdof)-2] .= newdof[1:end-1]
            else
                dof[t][j:j+length(newdof)-2] .= reverse(newdof[2:end])
            end
            j += length(newdof) - 1
        end
        dof[t][j:end] = k:k+(length(dof[t])-j)
        k += length(dof[t]) - j + 1
    end
    return dof
end


"""
    marked_dof(mesh::MeshHP{F,I,P},marker) where {F,I,P}

Returs a list of indices corresponding to the degrees of freedom marked with `marker`. If a vector of markers is passed, it returs degrees of freedom marked with any of them. 

Internally, `marked_dof` 
"""
function marked_dof(mesh::MeshHP{F,I,P}, marker) where {F,I,P}
    by_edge = degrees_of_freedom_by_edge(mesh)
    marked_dof(mesh, by_edge, marker)
end
function marked_dof(mesh::MeshHP{F,I,P}, by_edge, marker::N) where {F,I,P,N<:Integer}
    marked_dof(mesh, by_edge, [marker])
end
function marked_dof(mesh::MeshHP{F,I,P}, by_edge, markerslist::AbstractVector) where {F,I,P}
    (; edgelist) = mesh
    v = I[]
    for e in edges(edgelist)
        if marker(edgelist[e]) in markerslist
            push!(v, by_edge[e]...)
        end
    end
    return unique(v)
end

"""
    boundary_dof(mesh::MeshHP{F,I,P}) where {F,I,P}

Returns a list of the degrees of freedom lying at the boundary of the mesh.
Internally, it calls `degrees_of_freedom_by_edge` in order to obtain the nodal degrees of freedom. If a dof by edge dictionary has already been computed, it is recommended to run: 

    boundary_dof(mesh,by_edge) 
"""
function boudary_dof(mesh::MeshHP{F,I,P}) where {F,I,P}
    by_edge = degrees_of_freedom_by_edge(mesh)
    boundary_dof(mesh, by_edge)
end
function boundary_dof(mesh::MeshHP{F,I,P}, by_edge) where {F,I,P}
    marked_dof(mesh, by_edge, [1, 2])
end

function dirichlet_dof(mesh::MeshHP{F,I,P}) where {F,I,P}
    marked_dof(mesh, 1)
end



function transform_matrix!(A, vert)
    @views A[:, 1] .= 0.5(vert[:, 3] - vert[:, 2])
    @views A[:, 2] .= 0.5(vert[:, 1] - vert[:, 2])
end

function transform_term!(b, vert)
    @views b .= 0.5(vert[:, 1] + vert[:, 3])
end

"""
    hat_mass(B::Basis)

computes the mass matrix in the reference triangle fot the basis `B`. 
"""
function hat_mass(B::Basis)
    (; dim, b, C) = B
    sch = grundmann_moeller(Float64, Val(2), 2dim + 1)
    n = length(b)
    M = zeros(n, n)
    T = [[-1, -1], [1, -1], [-1, 1]]
    @inbounds for j = 1:n, i = 1:j
        M[i, j] = integrate(x -> b[i](x) * b[j](x), sch, T)
    end
    return SMatrix{n,n}(C' * Symmetric(M) * C)
end;



"""
    hat_stiff(B::Basis)

computes the stiff matrix in the reference triangle fot the basis `B`. 
"""
function hat_stiff(base::Basis)
    (; dim, ∇b) = base
    sch = grundmann_moeller(Float64, Val(2), 4 * dim + 1)
    S = zeros(dim, dim, 4)
    t̂ = [[-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]
    @inbounds for j = 1:dim, i = 1:j
        S[i, j, 1] = integrate(x -> ∇b[i](x)[1] * ∇b[j](x)[1], sch, t̂)
        S[i, j, 2] = integrate(x -> ∇b[i](x)[2] * ∇b[j](x)[1], sch, t̂)
        S[i, j, 3] = integrate(x -> ∇b[i](x)[1] * ∇b[j](x)[2], sch, t̂)
        S[i, j, 4] = integrate(x -> ∇b[i](x)[2] * ∇b[j](x)[2], sch, t̂)
    end
    @inbounds for k = 1:4
        S[:, :, k] = Symmetric(S[:, :, k])
    end
    SArray{Tuple{dim,dim,4},Float64}(S)
end;

"""
    hat_conv(base::Basis)
computes the matrix corresponding  
"""
function hat_conv(base::Basis)
    (; dim, b, ∇b) = base
    sch = grundmann_moeller(Float64, Val(2), 4 * dim + 1)
    C₁ = zeros(dim, dim)
    C₂ = zeros(dim, dim)
    t̂ = [[-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]
    @inbounds for j = 1:dim, i = 1:dim
        C₁[i, j] = integrate(x -> ∇b[i](x)[1] * b[j](x), sch, t̂)
        C₂[i, j] = integrate(x -> ∇b[i](x)[2] * b[j](x), sch, t̂)
    end
    (SMatrix{(dim, dim),Float64}(C₁), SMatrix{(dim, dim),Float64}(C₂))
end


"""
mass
Computes the mass matrix for a given mesh and a set of DOFs.

Parameters
----------
mesh : MeshHP
    The mesh object.
dof : Array{Int, 2}
    The set of DOFs
"""
function mass(mesh::MeshHP, dof, MD::MatrixDict)
    (; points, trilist) = mesh
    ℓ = maximum(maximum.(dof))
    I = Vector{Int64}(undef, ℓ)
    J = Vector{Int64}(undef, ℓ)
    V = Vector{Float64}(undef, ℓ)
    Aₜ = MMatrix{2,2}(zeros(2, 2))
    r = 1
    @inbounds for t in trilist
        dofT = dof[t]
        ℓₜ = length(dofT)
        p, pnod = pnodes(t, mesh)
        transform_matrix!(Aₜ, view(points, :, pnod))
        dAₜ = abs(det(Aₜ))
        n = compute_dimension(p)
        v = dAₜ * MD[p]
        i = repeat(dofT, n)
        j = repeat(dofT, inner = n)
        I[r:r+ℓₜ^2-1] = i
        J[r:r+ℓₜ^2-1] = j
        V[r:r+ℓₜ^2-1] = v
        r += ℓₜ^2
    end
    sparse(I, J, V)
end


"""
    stiff(mesh::MeshHO,dof,sd::ArrayDict,bd::BasisDict)

computes the stiffness matrix for a given `mesh`. It uses the integrals of every local basis present in the mesh corresponding to the reference triangle and stored in `sd`. 
"""
function stiff(mesh::MeshHP, dof, sd::ArrayDict, bd::BasisDict)
    (; points, trilist) = mesh
    ℓ = sum(x -> length(x)^2, dof)
    #ℓ = maximum(maximum.(dof))
    I = Vector{Int64}(undef, ℓ)
    J = Vector{Int64}(undef, ℓ)
    V = Vector{Float64}(undef, ℓ)
    Aₜ = MMatrix{2,2}(zeros(2, 2))
    iAₜ = MMatrix{2,2}(zeros(2, 2))
    # iAₜ² = MMatrix{2,2}(zeros(2,2))
    r = 1
    @inbounds for t in triangles(trilist)
        dofT = dof[t]
        p, pnod = pnodes(t, mesh)
        (; dim, C) = bd[p]
        transform_matrix!(Aₜ, view(points, :, pnod))
        iAₜ .= inv(Aₜ)
        iAₜ .= iAₜ * iAₜ'
        z = vec(iAₜ)
        dAₜ = abs(det(Aₜ))
        v = zeros(dim, dim)
        S = sd[p]
        for j = 1:dim, i = 1:j
            v[i, j] = S[i, j, :] ⋅ z
        end
        v = dAₜ * C' * Symmetric(v) * C
        i = repeat(dofT, dim)
        j = repeat(dofT, inner = dim)
        I[r:r+dim^2-1] = i
        J[r:r+dim^2-1] = j
        V[r:r+dim^2-1] = v
        r += dim^2
    end
    sparse(I, J, V)
end

"""
    conv(v::AbstractVector,mesh::MeshHO,dof,cd::BiMatrixDict,bd::BasisDict)

computes the convection matrix for a given velocity `v⃗` over a given `mesh`. It uses the integrals of every local basis present in the mesh corresponding to the reference triangle and stored in `cd`. 
"""
function conv(v::AbstractVector, mesh::MeshHP, dof, cd::BiMatrixDict, bd::BasisDict)
    (; points, trilist) = mesh
    ℓ = sum(x -> length(x)^2, dof)
    #ℓ = maximum(maximum.(dof))
    I = Vector{Int64}(undef, ℓ)
    J = Vector{Int64}(undef, ℓ)
    V = Vector{Float64}(undef, ℓ)
    Aₜ = MMatrix{2,2}(zeros(2, 2))
    iAₜ = MMatrix{2,2}(zeros(2, 2))
    # iAₜ² = MMatrix{2,2}(zeros(2,2))
    r = 1
    @inbounds for t in triangles(trilist)
        dofT = dof[t]
        p, pnod = pnodes(t, mesh)
        (; dim, C) = bd[p]
        transform_matrix!(Aₜ, view(points, :, pnod))
        iAₜ .= inv(Aₜ)'
        dAₜ = abs(det(Aₜ))
        C₁, C₂ = cd[p]
        z = zeros(dim, dim)
        @inbounds for j = 1:dim, i = 1:j
            z[i, j] = v' * iAₜ' * [view(C₁, i, j), view(C₂, i, j)]
        end
        v = dAₜ * C' * z * C
        I[r:r+dim^2-1] = repeat(dofT, dim)
        J[r:r+dim^2-1] = repeat(dofT, inner = dim)
        V[r:r+dim^2-1] = vec(v)
        r += dim^2
    end
    sparse(I, J, V)
end

""" 
    rhs(mesh,BD,f)

computes the right-hand side, ∫ϕf corresponing to `mesh` with data `f`. 
"""
function rhs(mesh, dof, BD, f)
    (; points, edgelist) = mesh
    ℓ = sum(length, dof)
    pmax = maximum(degree.(edgelist))
    #Constructors
    J = Vector{Int64}(undef, ℓ)
    V = Vector{Float64}(undef, ℓ)
    Aₜ = MMatrix{2,2}(zeros(2, 2))
    bₜ = MVector{2}(zeros(2))
    t̂ = [[-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]
    sch = grundmann_moeller(Float64, Val(2), 4pmax + 1)
    r = 1
    @inbounds for t in triangles(mesh)
        p, pnod = pnodes(t, mesh)
        transform_matrix!(Aₜ, view(points, :, pnod))
        transform_term!(bₜ, view(points, :, pnod))
        dAₜ = abs(det(Aₜ))
        (; dim, b, C) = BD[p]
        F = zeros(dim)
        for i in eachindex(F)
            F[i] = dAₜ * integrate(x -> b[i](x) * f(Aₜ * x + bₜ), sch, t̂)
        end
        v = C' * F
        J[r:r+dim-1] .= dof[t]
        V[r:r+dim-1] .= v
        r += dim
    end
    Vector(sparsevec(I, V, maximum(maximum.(dof))))
end




#Dictionaries

"""
    basisdict(mesh::MeshHP{F,I,P}) where {F,I,P}

constructs the `BasisDict` for a given `mesh`.
"""
function basisdict(mesh::MeshHP{F,I,P}) where {F,I,P}
    bd = BasisDict{P}()
    update!(bd, mesh)
end

"""
    update!(bd::BasisDict{P},mesh::MeshHP{F,I,P}) where {F,I,P}

updates the `BasisDict` `bd` for the given `mesh`, by adding all new `DegTuples`.
"""
function update!(bd::BasisDict{P}, mesh::MeshHP{F,I,P}) where {F,I,P}
    for t in triangles(mesh)
        p, _ = pnodes(t, mesh)
        isin, _ = gettoken(bd, p)
        if !isin
            set!(bd, p, Basis(p))
        end
    end
    return bd
end


"""
    stiffdict(mesh::MeshHP{F,I,P},bd::BasisDict) where {F,I,P}

constructs an `ArrayDict` storing the integrals of the form `∫∂Φᵢ∂Φⱼ` over the reference triangle.
"""
function stiffdict(bd::BasisDict{P}) where {P}
    sd = ArrayDict{P}()
    update!(sd, bd)
end
function update!(sd::ArrayDict, bd::BasisDict)
    @inbounds for p in keys(bd)
        isin, _ = gettoken(sd, p)
        if !isin
            set!(sd, p, hat_stiff(bd[p]))
        end
    end
    return sd
end


"""
    massdict(bd::BasisDict{P}) where {P}

constructs a dictionary with the same keys as `bd` but 
"""
function massdict(bd::BasisDict{P}) where {P}
    md = MatrixDict{P}()
    update!(md, bd)
end

function update!(md::MatrixDict{P}, bd::BasisDict{P}) where {P}
    @inbounds for p in keys(bd)
        isin, _ = gettoken(md, p)
        if !isin
            set!(md, p, hat_mass(bd[p]))
        end
    end
    return md
end

function convdict(bd::BasisDict{P}) where {P}
    cd = BiMatrixDict{P}()
    update!(cd, bd)
end

function update!(cd::BiMatrixDict{P}, bd::BasisDict{P}) where {P}
    @inbounds for p in keys(bd)
        isin, _ = gettoken(cd, p)
        if !isin
            set!(cd, p, hat_conv(bd[p]))
        end
    end
    return cd
end
