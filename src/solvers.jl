include("matrices.jl")
using LinearSolve
using CuthillMcKee


function poisson_rectangle(a,b,h,f)
    mesh   = rectangular_mesh(a,b,h)
    bd     = BasisDict()
    compute_pT!(mesh,bd)
    compute_dof!(mesh,bd)
    SD     = MatrixDict()
    for p in mesh.pE
        SD[p] = compute_hat_stiff(bd[p])
    end
    S      = compute_stiff(mesh,SD)
end

function laplacian_rectangle(a,b,h,p,f;method=:inexact)
    mesh       = rectangular_mesh(a,b,h)
    dof,dofT,ndof   = dof_by_triangle_eq(mesh,p)
    idof       = dof_interior(dof,a,b)
    S          = global_stiff_eq(mesh,p,dofT)[idof,idof]
    b          = global_rhs_eq(mesh,p,dof,dofT,f)[idof]
    u          = zeros(length(dof))
    perm       = 1:length(b)#symrcm(S)
    iperm      = 1:length(b)#symrcm(S,true,true)
    linprob    = LinearProblem(S[perm,perm],Array(b)[perm],alias_A=true,alias_b=true)
    init(linprob)
    sol        = solve(linprob)
    u[idof]   .= sol.u[iperm]
    
    #u[idof]    .= cg(S[idof,idof],b[idof],rtol=1e-15,atol=1e-15)[1]
    return mesh,dofT,dof,ndof,u
    #show_sol(dof,ndof,u)
end;

function error_L²(mesh::HPMesh,bd::BasisDict,u::Vector{Float64},sol::Function)
    (;P,T,pT,pE,Da) = mesh
    pmax  = maximum(pE)
    sch   = grundmann_moeller(Float64,Val(2),4pmax+1)
    t̂     = [[-1.,-1.],[1.,-1.],[-1.,1.]]
    err   = 0.
    Aₜ    = MMatrix{2,2}(zeros(2,2))
    bₜ    = MVector{2}(zeros(2))
    @inbounds for k in 1:size(T,2)
        p      = NTuple{3,Int8}(pT[:,k])
        t      = view(T,:,k)
        Aₜ    .= transform_matrix(P,t) 
        bₜ    .= transform_term(P,t)
        dAₜ    = abs(det(Aₜ))
        B      = bd[p]
        (;b,C) = B
        U(x)   = u[Da[k]]⋅(C'*[φ(x) for φ in b])
        err += dAₜ*integrate(x->(U(x)-sol(Aₜ*x+bₜ))^2,sch,t̂)  
    end
    √err
end

function error_L²_eq(u,dof,dof_T,p,sol)
    #sch   = grundmann_moeller(Float64,Val(2),4p+1)
    B     = standard_basis(p)
    T̂     = SMatrix{3,2,Float64}([-1 -1;1 -1;-1 1])
    X,W   = simplexquad(4p,T̂)
    X    = [X[i,:] for i in 1:size(X,1)]
    nodes = boundary_nodes(p)
    C     = matrix_C(B,nodes)
    error = 0.
    Aₜ     = MMatrix{2,2}(zeros(2,2))
    bₜ     = MVector{2}(zeros(2))
    for t in dof_T
        T    = dof[[t[1],t[p+1],t[2p+1]]]
        Aₜ   .= 0.5*[T[2]-T[1] T[3]-T[1]] 
        bₜ   .= (T[2]+T[3])/2
        dAₜ   = abs(det(Aₜ))
        U(x) = u[t]⋅(C'*[B[j](x) for j in 1:length(t)])
        #error += dAₜ*integrate(x->(U(x)-sol(Aₜ*x+bₜ))^2,sch,T̂')
        error += dAₜ*integral(x->(U(x)-sol(Aₜ*x+bₜ))^2,X,W)
    end
    return √error
end;


function eigenvalues_rectangle(a,b,h₁,h₂,p,nv)
    (mesh,vor)    = rectangle_mesh(a,b,h₁,h₂)
    dof,dofT,ndof = dof_by_triangle_eq(mesh,p)
    idof          = .!boundary_dof_rectangle(a,b,dof)
    S             = global_stiff_eq(mesh,p,dofT)
    M             = global_mass_eq(mesh,p,dofT)
    λ,x,info      = geneigsolve((-S[idof,idof],-M[idof,idof]),howmany=nv,which=:SR,issymetric=true,isposdef=true)
    u             = zeros(length(dof),nv)
    u[idof,:]     .= x
    return mesh,dof,ndof,λ,u
end;

function show_sol(dof,ndof,z)
    X = hcat(dof[ndof]...)
    return surface(X[1,:],X[2,:],z[ndof])
end;
