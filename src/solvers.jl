"""
    DegreesOfFreedom(mesh::MeshHP)

returns a `DegreesOfFreedom{I}` structure that stores the degrees of freedom of a given mesh. Two dictionaries (from `Dictionaries.jl`) are used: 
    + `by_edge`: the keys are the edges of the mesh and the values are vectors of indices for the degrees of freedom. This dictionary contains all the nodal dof.
    + `by_tri`: similar to `by_edge`, but including the indices of the degrees of freedom corresponding to bubble functions. 
The type `I<:Integer` is the type used for indexing triangles and edges in `mesh`.
"""
mutable struct DegreesOfFreedom{I<:Integer}
    by_tri::DofDict{I}
    by_edge::EdgeDofDict{I}
end

function DegreesOfFreedom(mesh::MeshHP)
    by_edge = degrees_of_freedom_by_edge(mesh)
    by_tri = degrees_of_freedom(mesh, by_edge)
    DegreesOfFreedom(by_tri, by_edge)
end

"""
    BoundaryConditions()

Returns a `BoundaryConditions` structure that stores the functions corresponding to the `dirichlet` and `neumann` conditions.

    BoundaryConditions(f::Function)
sets `dirichlet = f` and `neumann = nothing`.

    Named arguments can be used:

```julia

    julia> f(x) = x[1] > 0 ? 1. : 0.;
    julia> g(x) = x[1]^2 + x[2]^2
    julia> bc = BoundaryConditions(dirichlet = f, neumann = g)
```
"""
struct BoundaryConditions
    dirichlet::Union{Function,Nothing}
    neumann::Union{Function,Nothing}
end
BoundaryConditions(f::Function) = BoundaryConditions(f, nothing)
BoundaryConditions(; dirichlet = isequal(0.0), neumann = isequal(0.0)) =
    BoundaryConditions(dirichlet, neumann)

"""
    ReferenceDicts(mesh)

Returns a `ReferenceDicts` structure that stores dictionaries (from `Dictionaries.jl`) containing information computed on the reference triangle `T̂`, with vertices `[-1,1]`, `[-1,-1]` and `[1,-1]`. All dictionaries have keys given by `DegTuple`s corresponding to the combinations of degrees `(p₁,p₂,p₃)` that are present in a given mesh.  `P` is the type of the degrees. 
The dictionaries are: 
    + `basisref`: stores the al the information corresponding to the local basis and its gradientes.
    + `stiffref`: stores integrals of the form `∫∇Φᵢ∇ΦⱼdT̂` that are used for the assembly of the stiff matrix. 
    + `massref`: stores integrals of the form `∫ΦᵢΦⱼdT̂` that are used for the assembly of the mass matrix.
    + `convref`: stores integrals of the form `∫∇ΦᵢΦⱼdT̂` that are used for the assembly of the convection matrix. 
If some matrix is not present on an instance of `ConstantCoeffProblem`, the corresponding dictionary is set to `nothing`. 
"""
struct ReferenceDicts{P<:Integer}
    basisref::Union{BasisDict{P},Nothing}
    stiffref::Union{ArrayDict{P},Nothing}
    massref::Union{MatrixDict{P},Nothing}
    convref::Union{BiMatrixDict{P},Nothing}
end
function ReferenceDicts(mesh::MeshHP{F,I,P}) where {F,I,P}
    bd = basisdict(mesh)
    sd = stiffdict(bd)
    md = massdict(mesh, bd)
    cd = convdict(mesh, bd)
    ReferenceDicts(bd, sd, md, cd)
end


"""
    update!(ref::ReferenceDicts,mesh)

updates the dictionaries stored in `ref` according to the new information present in `mesh`. In particular, if a new `DegTuple` appear in `mesh` after `p`-refinement, `update!` will compute all the relevant information for computing and assembling matrices accordingly. 
"""
function update!(ref::ReferenceDicts, mesh)
    (; bd, sd, md, cd) = ref
    update!(bd, mesh)
    update!(sd, mesh)
    if md !== nothing
        update!(md, mesh)
    end
    if cd !== nothing
        update!(cd, mesh)
    end
end

"""
    ConstantCoeff{F<:AbstractFloat,U<:Function}

Stores the coefficients and right hand side of a `ConstantCoeffProblem`. In other words, for a problem of the form:  

    α∫∇u⋅∇Φ + ∫v⃗⋅∇u Φ + c∫uΦ = ∫fΦ
the `ConstantCoeff` `struct` stores the numbers `α` and `c`, the vector `v⃗` and the function `f`. 

    ConstantCoeff(α,v⃗,c,f)
creates a `struct` of constant coefficients promoting `α,c` and the coeffients of `v⃗` to a common number type. 

    ConstantCoeff{F}(α,v⃗,c,f) where F<:AbstractFloat
creates a `struct` of constant coefficients converting all coefficients to type `F`.

    ConstantCoeff(α,c,f)
sets `v⃗` to a vector of zeros. 

    ConstantCoeff(α,f)
sets `v⃗` and `c` to zero.

    ConstantCoeff(α)
sets `v⃗` and `c` to zero and `f=isequal(0)`. 
"""
struct ConstantCoeff{F<:AbstractFloat,U<:Function}
    α::F
    v⃗::SVector{2,F}
    c::F
    f::U
end
function ConstantCoeff(α, v⃗, c, f)
    α, v₁, v₂, c = promote(α, v⃗..., c)
    v⃗ = SVector{2}([v₁, v₂])
    ConstantCoeff(α, v⃗, c, f)
end
ConstantCoeff{F}(α, v⃗, c, f) where {F<:AbstractFloat} = ConstantCoeff(F(α), SVector{2,F}(v⃗), F(c), f)
ConstantCoeff(α, v⃗::V, f) where {V<:AbstractVector} = ConstantCoeff(α, v⃗, convert(typeof(α), 0), f)
ConstantCoeff(α, c, f) = ConstantCoeff(α, SVector{2,typeof(α)}(zeros(2)), c, f)
ConstantCoeff(α, f::U) where {U<:Function} = ConstantCoeff(α, convert(typeof(α), 0), f)
ConstantCoeff(α) = ConstantCoeff(α, isequal(0))

"""
   ConstantCoeffProblem(mesh,α,v⃗,c,f,bc)

creates problem of the form 
    α∫∇u⋅∇ϕ + ∫v⃗⋅∇u ϕ + c∫uϕ = ∫fϕ
where `mesh` is the mesh, `α`, `v⃗` and `c` are constants, `f` is a given function and `bc` are boundary conditions.  
"""
mutable struct ConstantCoeffProblem{F<:AbstractFloat,I<:Integer,P<:Integer}
    Ω::MeshHP{F,I,P}
    coeff::ConstantCoeff{F}
    ref::ReferenceDicts{P}
    dof::Union{DegreesOfFreedom{I},Nothing}
    bc::BoundaryConditions
end

#General constructor
function ConstantCoeffProblem(
    Ω::MeshHP{F,I,P},
    α,
    v⃗::AbstractVector,
    c,
    f::Function,
    bc::BoundaryConditions,
) where {F,I,P}
    bd = basisdict(Ω)
    sd = stiffdict(bd)
    md = c != 0 ? massdict(Ω, bd) : nothing
    cd = v⃗[1] != 0 || v⃗[2] != 0 ? convdict(Ω, bd) : nothing
    ref = ReferenceDicts(bd, sd, md, cd)
    dof = DegreesOfFreedom(Ω)
    coeff = ConstantCoeff(α, SVector{2,F}(v⃗), c, convert(F, c), f)
    ConstantCoeffProblem(Ω, coeff, ref, dof, bc)
end

#Constructor for diffusion
function ConstantCoeffProblem(
    Ω::MeshHP{F,I,P},
    α::R,
    f::Function,
    bc::BoundaryConditions,
) where {F,I,P,R<:Number}
    ConstantCoeffProblem(Ω, F(α), SVector{2,F}(zeros(2)), F(0), f, bc)
end
#Constructor for reaction-diffusion
function ConstantCoeffProblem(
    Ω::MeshHP{F,I,P},
    α::R,
    c::R,
    f::Function,
    bc::BoundaryConditions,
) where {F,I,P,R<:Number}
    ConstantCoeffProblem(Ω, F(α), SVector{2,F}(zeros(2)), F(c), f, bc)
end

function ConstantCoeffProblem(
    Ω::MeshHP{F,I,P},
    α::R,
    v⃗::V,
    f::Function,
    bc::BoundaryConditions,
) where {F,I,P,R<:Number,V<:AbstractVector}
    ConstantCoeffProblem(Ω, F(α), SVector{2,F}(v⃗), F(0), f, bc)
end

function update!(prob::ConstantCoeffProblem)
    (; Ω, ref, dof) = prob
    dof = DegreesOfFreedom(Ω)
    update!(ref)
end



function CommonSolve.solve(prob::ConstantCoeffProblem)
    (; Ω, coeff, ref, dof, bc) = prob
    (; α, v⃗, c, f) = coeff
    (; basisref, stiffref, massref, convref) = ref
    (; by_tri, by_edge) = dof
    A = stiff(Ω, by_tri, stiffref, basisref)
    if c ≠ 0
        A += mass(Ω, by_tri, massref, basisref)
    end
    if v⃗[1] ≠ 0 || v⃗[1] ≠ 0
        A += conv(Ω, by_tri, convref, basisref)
    end
    b = rhs(Ω, dof, bd, f)
    u = zeros(length(b))
    ∂dof = boundary_dof(Ω, by_edge)
    idof = setdiff(1:maximum(maximum.(by_tri)), ∂dof)
    linprob = LinearProblem(A[idof, idof], b[idof])#,alias_A=true,alias_b=true)
    sol = solve(linprob, abstol = 1e-12, reltol = 1e-14)
    u[idof] = sol.u #A[idof,idof]\b[idof]
    return u
end



# function error_L²(prob::ConstantCoeffProblem,u::Vector{Float64},sol::Function)
#     (;Ω,bd) = prob
#     (;points,trilist,edgelist) = Ω
#     pmax  = maximum(degree.(edgelist))
#     sch   = grundmann_moeller(Float64,Val(2),4pmax+1)
#     t̂     = [[-1.,1.],[-1.,-1.],[1.,-1.]]
#     err   = 0.
#     Aₜ    = MMatrix{2,2}(zeros(2,2))
#     bₜ    = MVector{2}(zeros(2))
#     dof   = degrees_of_freedom(Ω)
#     @inbounds for t in triangles(trilist)
#         p,pnod    = pnodes(t,Ω)
#         transform_matrix!(Aₜ,view(points,:,pnod)) 
#         transform_term!(bₜ,view(points,:,pnod))
#         dAₜ    = abs(det(Aₜ))
#         (;b,C) = bd[p]
#         U(x)   = u[dof[t]]⋅(C'*[φ(x) for φ in b])
#         err += dAₜ*integrate(x->(U(x)-sol(Aₜ*x+bₜ))^2,sch,t̂)  
#     end
#     √err
# end

# function error_L²_eq(u,dof,dof_T,p,sol)
#     #sch   = grundmann_moeller(Float64,Val(2),4p+1)
#     B     = standard_basis(p)
#     T̂     = SMatrix{3,2,Float64}([-1 -1;1 -1;-1 1])
#     X,W   = simplexquad(4p,T̂)
#     X    = [X[i,:] for i in 1:size(X,1)]
#     nodes = boundary_nodes(p)
#     C     = matrix_C(B,nodes)
#     error = 0.
#     Aₜ     = MMatrix{2,2}(zeros(2,2))
#     bₜ     = MVector{2}(zeros(2))
#     for t in dof_T
#         T    = dof[[t[1],t[p+1],t[2p+1]]]
#         Aₜ   .= 0.5*[T[2]-T[1] T[3]-T[1]] 
#         bₜ   .= (T[2]+T[3])/2
#         dAₜ   = abs(det(Aₜ))
#         U(x) = u[t]⋅(C'*[B[j](x) for j in 1:length(t)])
#         #error += dAₜ*integrate(x->(U(x)-sol(Aₜ*x+bₜ))^2,sch,T̂')
#         error += dAₜ*integral(x->(U(x)-sol(Aₜ*x+bₜ))^2,X,W)
#     end
#     return √error
# end;


# function eigenvalues_rectangle(a,b,h₁,h₂,p,nv)
#     (mesh,vor)    = rectangle_mesh(a,b,h₁,h₂)
#     dof,dofT,ndof = dof_by_triangle_eq(mesh,p)
#     idof          = .!boundary_dof_rectangle(a,b,dof)
#     S             = global_stiff_eq(mesh,p,dofT)
#     M             = global_mass_eq(mesh,p,dofT)
#     λ,x,info      = geneigsolve((-S[idof,idof],-M[idof,idof]),howmany=nv,which=:SR,issymetric=true,isposdef=true)
#     u             = zeros(length(dof),nv)
#     u[idof,:]     .= x
#     return mesh,dof,ndof,λ,u
# end;

# function show_sol(dof,ndof,z)
#     X = hcat(dof[ndof]...)
#     return surface(X[1,:],X[2,:],z[ndof])
# end;
