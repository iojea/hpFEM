include("basis.jl")
using Base.Iterators
using Plots
import Plots.plot
using ColorSchemes
using Triangulate
using StatsBase
using UnPack
using Printf





"""

  HPMesh

is the struct that stores the mesh. Fields are: 

+ `P::Matrix{<:AbstractFloat}`: matrix of points, size: `2×N`.
+ `T::Matrix{U1} where U1<:Unsigned`: matrix of triangles, size `3×NT`. The indices in each column are sorted according the degree attached to the edge, i.e.: `(T[1,k],T[2,k])` is the edge with the smaller degree. 
+ `E::Matrix{U1} where U1<:Unsigned`: matrix of edges
+ `ET::Matrix{U1} where U1<:Unsigned`: matrix with the edges of each triangle, size: `3xNT`
+ `Pm::Vector{U2} where U2<:Unsigned`: point markers
+ `Em::Vector{U2} where U2<:Unsigned`: edge markers
+ `pE::Vector{U2} where U2<:Unsigned`: value of `p` on each edge
+ `pT::Matrix{U2} where U2<:Unsigned`: matrix of values of `p` on each triangle (sorted). Size `3×NT`
+ `D::Matrix{<:AbstractFloat}`: matrix of degrees of freedom, size: `2×ND`. Degrees of freedom corresponding to bubble functions are placed at the barycenter of the triangle. 
+ `Da::Vector{Vector{U1} where U1<:Unsigned}`: list of degrees of freedon on each triangle.
+ `rT::Matrix{U2} whew U2<:Unsigned`: edges of each triangle marked for refinement. The mark uses Verfürth's code: 
    - `0`: no refinement.
    - `i`: refine edge i (1<=i<=3)
    - `4`: red refinement
    - `4k+i`: blue refinement of edges `i` and `i+k`
+ `rE::Vector{Bool}`: edges marked for refinement.
+ `leT::Vector{U2} where U2<:Unsigned`: longest edge of triangle `T`. 

Typically, `U1=UInt32`, `U2=UInt8`, `AbstractFloat->Float64`.
"""
struct HPMesh{#basic types
              F<:AbstractFloat, #Float64
              I<:Integer, #UInt32
              U<:Integer, #UInt8
              #types of variables
              MF<:AbstractMatrix{F},
              MI<:AbstractMatrix{I},
              VU<:AbstractVector{U},
              MU<:AbstractMatrix{U},
              VB<:Vector{Bool},
              VVR<:Union{Nothing,AbstractVector{VF} where VF<:AbstractVector{F}},
              VVI<:Union{Nothing,AbstractVector{VI} where VI<:AbstractVector{I}}
             }
    P::MF
    T::MI
    E::MI
    ET::MI
    Pm::VU
    Em::VU
    pE::VU
    pT::MU
    rT::VU
    rE::VB
    leT::VU
    D::VVR
    Da::VVI
end;
#HPMesh(P::Matrix{F},T::Matrix{I1},E::Matrix{I1},ET::Matrix{I1},Pm::Vector{I2},Em::Vector{I2},pE::Vector{I2},pT::Matrix{I2},rT::Vector{I2},rE::Vector{Bool},leT::Vector{I2}) where {F<:AbstractFloat,I1<:Integer,I2<:Integer} = HPMesh(P,T,E,ET,Pm,Em,pE,pT,rT,rE,leT,nothing,nothing)

"""

    circular_mesh(c,r,h)

builds a circular mesh with center `c`, radious `r` and edge-length `h`.
"""
function circular_mesh(c,r,h)
    tri       = TriangulateIO()
    n         = Int(1+2π*r÷h)
    θ         = range(start=0,stop=2π,length=n)[1:end-1]
    pointlist = hcat([[c[1]+cos(t),c[2]+sin(t)] for t in θ]...)
    edgelist  = hcat([[i,i+1] for i in 1:length(pointlist)-1]...,
                  [length(pointlist),1])
    edgemarkerlist = ones(length(edgelist))
    @pack! tri = pointlist,edgelist,edgemarkerlist
    maxa = Printf.@sprintf "%.15f" h^2/2
    angle= Printf.@sprintf "%.15f" 30.
    (tri,vor) = triangulate("ea$(maxa)q$(angle)Qv",tri)

    args  = [tri.pointlist, #P
             UInt32.(tri.trianglelist), #T
             UInt32.(tri.edgelist),#E
             compute_ET(size(tri.trianglelist,2),vor.edgelist),#ET
             UInt8.(tri.pointmarkerlist),#Pm#
             UInt8.(tri.edgemarkerlist),#Em
             ones(UInt8,length(tri.edgemarkerlist)),#pE
             ones(UInt8,size(tri.trianglelist)...),#pT
             zeros(UInt8,size(tri.trianglelist,2)),#rT
             zeros(Bool,length(tri.edgemarkerlist)),#rE
             ones(UInt8,size(tri.trianglelist,2)),#leT
             nothing,
             nothing]
    println.(typeof.(args))
    mesh = HPMesh(args...)
    compute_leT!(mesh)
    mesh
end;
circular_mesh(h::T) where T<:AbstractFloat = circular_mesh(zeros(2),1.,h) 

"""
    rectangular_mesh(a,b,h)

builds a mesh of the rectangle `[0,a]×[0,b]`, with edge-length `h`. The values of `p` are initialized as `1`.
"""
function rectangular_mesh(a,b,h)
    tri       = TriangulateIO()
    pointlist = hcat([[[a,y] for y in 0:h:b-h]...,
                      [[x,b] for x in a:-h:h]...,
                      [[0.,y] for y in b:-h:h]...,
                      [[x,0.] for x in 0.:h:a-h]...]...)
    edgelist  = hcat([[i,i+1] for i in 1:length(pointlist)-1]...,
                  [length(pointlist),1])
    edgemarkerlist = ones(length(edgelist))
    @pack! tri = pointlist,edgelist,edgemarkerlist
    maxa = Printf.@sprintf "%.15f" h^2
    (tri,vor) = triangulate("ea$(maxa)qQv",tri)
    mesh = HPMesh(tri.pointlist,
                  UInt32.(tri.trianglelist),
                  UInt32.(tri.edgelist),
                  compute_ET(size(tri.trianglelist,2),vor.edgelist),
                  UInt8.(tri.pointmarkerlist),
                  UInt8.(tri.edgemarkerlist),
                  ones(UInt8,length(tri.edgemarkerlist)),
                  ones(UInt8,size(tri.trianglelist)...),
                  zeros(UInt8,size(tri.trianglelist,2)),
                  zeros(Bool,length(tri.edgemarkerlist)),
                  ones(UInt8,size(tri.trianglelist,2)),
                  nothing,
                  nothing,
    )
end;
square_mesh(a,h) = rectangular_mesh(a,a,h);


@inline almost(x,y;tol=1e-14) = isapprox(x,y,atol=tol,rtol=0)
@inline function boundary_rectangle(dof,a,b;tol=1e-14)
    almost.(dof[1,:],0;tol=tol).||almost.(dof[1,:],a;tol=tol).||almost.(dof[2,:],0;tol=tol).||almost.(dof[2,:],b;tol=tol)
end


"""

    compute_ET(T,vE)

Computes the matrix `ET` of a mesh, given the matrix of its triangles `T` and the edges of the dual graph (Voronoi diagram), stored in `vE`.
"""
function compute_ET(ntri,vE)
    ET  = zeros(UInt32,3,ntri)
    for k in eachindex(eachcol(ET))
        ET[:,k] = findall(x->k in x,eachcol(vE))
    end
    return ET
end


"""
    find_sorted_edges(t,E)

Given a tringle t (three indices) and an edge matrix E, it returns the indices of the edges in the same order as they are implicitly given by t: (t[1],t[2]),(t[2],t[3]),(t[3],t[1]).
"""
function find_sorted_edges(t,E)
    x  = zeros(Int32,3)
    @inbounds for i in 1:3
        x[i] = findfirst(z->issubset([t[i],t[mod1(i+1,3)]],z),eachcol(E))
    end 
    return x
end

"""
    function compute_pT!(mesh::HPMesh)

        Taking the values from `mesh.pE`, computes the matrix `mesh.pT`. It also rearranges the triangles so that the first edge corresponds to the minimum value of p. If the pair (p₁,p₂) has not been computed yet, it is added to BD. 
""" 
function compute_pT!(mesh::HPMesh,BD::BasisDict)
    (;E,pE,pT) = mesh
    ntri = size(mesh.T,2)
    q⃗    = MVector{3}(zeros(Int8,3))
    @inbounds for j in 1:ntri
        t       = view(mesh.T,:,j)
        et      = find_sorted_edges(t,E)
        q⃗       = pE[et]
        sort_tri!(t,q⃗)
        pT[:,j] = q⃗
        key     = NTuple{3,Int8}(pT[:,j])
        if key ∉ keys(BD)
            BD[key] = Basis(key)
        end 
    end
end;


"""

    sort_tri!(t,q⃗)

Given a triangle `t=[t₁,t₂,t₃]` and a vector `q⃗=[q₁,q₂,q₃]`, such that `q₁` correspond to the edge `(t₁,t₂)`, etc., `sort_tris` generates a permutation of `t` such that the corresponding values of `q` are sorted.
"""
function sort_tri!(t,q⃗)
    i = 1
    flag = !issorted(q⃗)
    while flag && i<3
        circshift!(t,1)
        circshift!(q⃗,1)
        i += 1
        flag = !issorted(q⃗)
    end
    if flag
        permute!(t,[1,3,2])
        permute!(q⃗,[1,3,2])
        while !issorted(q⃗) && i<6
            circshift!(t,1)
            circshift!(q⃗,1)
            i += 1
            flag = !issorted(q⃗)
        end
    end   
    return nothing 
end

"""

    is_in(v,V[,tol])

Decides if `v` is a column of `V` (with tolerance `tol`). Returns the index of the column of `nothing`. 
"""
function is_in(v,V;tol=1e-14)
    f(u) = norm(v-u)<tol
    findfirst(f,V)
end


"""

    transform_matrix(P,t) 

computes the matrix that transforms the reference triangle into the triangle nodes `p[:,t]`.
"""
@inline transform_matrix(P,t) = 0.5*[P[:,t[3]]-P[:,t[2]] P[:,t[1]]-P[:,t[2]]];


"""

    transform_term(P,t) 

computes the affine term of the transformation that transforms the reference triangle into the triangle nodes `p[:,t]`.
"""
@inline transform_term(P,t)   = (P[:,t[1]]+P[:,t[3]])/2


"""

    fill_bd!(pT::Matrix{Int8},bd::BasisDict)

fills the Basis dictionary using the values of `pT`.
"""
function fill_bd!(pT::Matrix{Int8},bd::BasisDict)
    for p in eachcol(pT)
        ptuple = NTuple{3,Int8}(p)
        if ptuple ∉ keys(bd)
            bd[ptuple] = Basis(ptuple)
        end
    end    
end

"""

    compute_dof!(mesh::HPMesh,bd::BasisDict)

computes the degrees of freedom of `mesh`, assuming that the values stored in `mesh.pT` are correct. It uses `bd` as a dictionary of already available bases (only uses `bd[key].nodes`) 
"""
    
function compute_dof!(mesh::HPMesh
,bd::BasisDict)
    (;P,T,E,Em,pT,pE,D,Da) = mesh
    nP   = size(P,2)
    ndof_by_t = [bd[NTuple{3,Int8}(p)].dim for p in eachcol(pT)]
    Da  .= zeros.(Int,ndof_by_t)
    nT   = sum(ndof_by_t) #sum of dofs by triangle
    Ei   = Em.==0
    nre  = sum(p for p in pE[Ei])-sum(Ei) #number of repeated nodal dofs that are not vertices
    nrv  = sum(collect(values(countmap(T))).-1)
    nD   = nT-nrv-nre
    resize!(D,nD)
    Aₜ   = MMatrix{2,2}(zeros(2,2))
    bₜ   = MVector{2}(zeros(2))
    bar  = MVector{2}(zeros(2))
    ncolsD    = size(D,2)
    k₁        = nP+1
    for j in eachindex(Da)
        p    = NTuple{3,Int8}(pT[:,j])
        t    = view(T,:,j) #indices of vertices such that p increases on the edges
        vind = [1,p[1]+1,p[1]+p[2]+1]
        Da[j][vind] .= t  
        if max(p...)>1
            nodes = bd[p].nodes
            Aₜ   .= transform_matrix(P,t)
            bₜ   .= transform_term(P,t)
            for k in Iterators.filter(∉(vind),1:size(nodes,2))
                z = Aₜ*nodes[:,k] + bₜ
                i = is_in(z,D[1:nP+1:min(k₁,ncolsD)])
                if i === nothing
                    D[k₁] = z
                    i       = k₁
                    k₁     += 1
                else
                    i = i+nP
                end
                Da[j][k] = i
            end
            bar .= sum(P[:,t],dims=2)/3
            for k in sum(p)+1:length(Da[j])
                D[k₁]  = bar
                Da[j][k] = k₁ 
                k₁ += 1
            end
        end
    end
end     


     
"""

    plot(mesh::HPMesh,color_scheme::Symbol)

plots `mesh`. Border edges are thicker. Color indicates the degree asigned to the edge.
"""
function plot(mesh::HPMesh
;color_scheme::Symbol=:blues,cb=true,kwargs...)
    @unpack P,E,Em,pE = mesh
    nc      = Int(maximum(pE))
    mc      = Int(minimum(pE))
    pal     = palette(color_scheme,nc-mc+2)#nc-mc>0 ? palette(color_scheme,nc-mc+1) : palette(color_scheme,1)
    plt     = plot(bg=:gray30,legend=:none,border=:none;kwargs...)
    l = @layout [a{0.95w} b]
    @inbounds for i in eachindex(pE)
        e = E[:,i]
        x = P[:,e]
        if Em[i]==1
            plot!(plt,x[1,:],x[2,:],c=pal[pE[i]-mc+1],linewidth=4)
        else
            plot!(plt,x[1,:],x[2,:],c=pal[pE[i]-mc+1],linewidth=1)
        end
    end
    if cb
        #cl2       = nc-mc > 0 ? nc-mc : 2
        color_bar = heatmap(rand(2,2),clims=(mc,nc+1),framestyle=:none,c=pal,cbar=true,lims=(-1,0))  
        return plot(plt,color_bar,layout=l)
    else
        return plt
    end
end;
#plot(mesh::HPMesh) = plot(mesh,:blues)

function plot_marks(mesh::HPMesh
;kwargs...)
    @unpack P,E,T,rE,rT = mesh
    plt = plot(bg=:gray30,legend=:none,border=:none;kwargs...)
    @inbounds for i in eachindex(rE)
        x = P[:,E[:,i]]
        plot!(x[1,:],x[2,:],c= rE[i] ? :black : :cornflowerblue)
    end
    x = []; y = []; col = []
    for i in eachindex(rT)
        if rT[i]>0
            b = sum(P[:,T[:,i]],dims=2)/3
            push!(x,b[1])
            push!(y,b[2])
            if 1≤rT[i]≤3
                push!(col,:green)
            elseif rT[i]==4
                push!(col,:red)
            else
                push!(col,:blue)
            end
        end
    end
    scatter!(x,y,c=col,ms=2)
end

function plot_marks_only(mesh::HPMesh
;kwargs...)
    @unpack P,T,E,ET,rE,rT = mesh
    plt = plot(bg=:gray30,legend=:none,border=:none;kwargs...)
    x = []; y = []; col = []
    for i in eachindex(rT)
        if any(rE[ET[:,i]])
            for j in ET[:,i]
                pts = P[:,E[:,j]]
                plot!(pts[1,:],pts[2,:],c= rE[j] ? :black : :cornflowerblue)
            end
        end
        if rT[i]>0
            b = sum(P[:,T[:,i]],dims=2)/3
            push!(x,b[1])
            push!(y,b[2])
            if 1≤rT[i]≤3
                push!(col,:green)
            elseif rT[i]==4
                push!(col,:red)
            else
                push!(col,:blue)
            end
        end
    end
    scatter!(x,y,c=col,ms=2)
end
