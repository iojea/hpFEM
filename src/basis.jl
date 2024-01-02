"""
    Basis():

Defines the elements of a standard basis: 
+ `dim` is the number of functions in the basis (`==length(b)`).
+ `nodes` are the nodes at the boundary of the reference triangle
+ `b` is the basis defined in terms of Legendre polynomials
+ `∇b` contains the gradients of the functions in `b`
+ `C` is the matrix that makes the change of variables to the mixed basis. 

The constructor `Basis(p₁[,p₂,p₃])` builds the basis corresponding to the edges degrees `p₁`, `p₂` and `p₃`. As usual, if `p₃` is not given it is assumed that `p₃=p₂`, and if `p₂` and `p₃` are not given, it is assumed that `p₂=p₃=p₁`. 
"""
struct Basis{N1,N2,I<:Integer,T<:AbstractFloat,F<:Function,V<:Vector{F},M1<:SMatrix{2,N1}{T}, M2<:SMatrix{N1,N2}{T}} 
    dim::I
    nodes::M1
    b::V 
    ∇b::V 
    C::M2
end

function Basis(p₁,p₂,p₃)
    println("Creating basis: ($p₁,$p₂,$p₃)")
    dim   = compute_dimension(p₁,p₂,p₃)
    println("Dimension = $dim")
    nodes = boundary_nodes(p₁,p₂,p₃)
    b     = standard_basis(p₁,p₂,p₃)
    ∇b    = ∇standard_basis(p₁,p₂,p₃)
    C     = matrix_C(b,nodes)
    Basis(dim,nodes,b,∇b,C)
end
Basis(p₁,p₂)   = Basis(p₁,p₂,p₂)
Basis(p₁::T) where T<:Integer  = Basis(p₁,p₁) 
Basis(t::T) where T<:Tuple = Basis(t...)
"""
    const BasisDict

Defines a dictionary with an `NTuple{3,Int}` containing the values `(p₁,p₂,p₃)` as key and a `Basis` as value. 
"""
const BasisDict = Dict{NTuple{3,UInt8},Basis}




"""

    compute_dimension(p₁,p₂,p₃)
    compute_dimension(p₁,p₂)   
    compute_dimension(p₁)
    compute_dimension(set::Set)

Computes the dimension of the space 𝒫p₁p₂p₃. 
"""
compute_dimension(p₁,p₂,p₃) = sum(min(p₂,p₃-j)+1 for j in 0:p₁);
compute_dimension(p₁,p₂)    = compute_dimension(p₁,p₂,p₂)
compute_dimension(p₁)       = compute_dimension(p₁,p₁)
compute_dimension(t::T) where T<:Tuple = compute_dimension(t...)



"""
    
    boundary_nodes(p₁[,p₂,p₃])

Defines the nodes on the boundary of the reference triangle, with vertices: [-1,-1],[1,-1],[-1,1]. The degrees of the polyomial restricted to the edges are assumed to satisfy p₁≤p₂≤p₃. `boundary_nodes(p₁,p₂)` takes p₃=p₂ and `boundary_nodes(p₁)` takes p₃=p₂=p₁. 
"""
function boundary_nodes(p₁,p₂,p₃)
    seg₁   = range(start=1,stop=-1,length=p₁+1)
    seg₂   = range(start=-1,stop=1,length=p₂+1)
    seg₃   = range(start=-1,stop=1,length=p₃+1)
    nodes  = SMatrix{2,Int(p₁+p₂+p₃)}(hcat(
                    stack([-1,y] for y in seg₁[1:end-1]), 
                    stack([x,-1] for x in seg₂[1:end-1]),
                    stack([-z,z]  for z in seg₃[1:end-1])
                    ))
end;

boundary_nodes(p₁,p₂) = boundary_nodes(p₁,p₂,p₂);
boundary_nodes(p₁::T) where T<:Integer = boundary_nodes(p₁,p₁,p₁);
boundary_nodes(t::Tuple) = boundary_nodes(t...)


"""

    standard_basis(p₁[,p₂])

Builds the standard basis on T̂, meaning: bₖ(x) = Pᵢ(x₁)Pⱼ(x₂), where Pᵢ is the i-th polynomial of Legendre. It returns a vector of functions. 
""" 
function standard_basis(p₁,p₂,p₃)
    n = compute_dimension(p₁,p₂,p₃)
    B = Vector{Function}(undef,n)
    k = 1
    for i in 0:p₁
        for j in 0:min(p₂,p₃-i)
            B[k] = x->Pl(x[1],j)*Pl(x[2],i)
            k   += 1
        end
    end
    return B
end;
standard_basis(p₁,p₂) = standard_basis(p₁,p₂,p₂)
standard_basis(p₁::T) where T<:Integer    = standard_basis(p₁,p₁);
standard_basis(t::Tuple) = standard_basis(t...)

"""

    ∇standard_basis(p₁[,p₂,p₃])

Builds the grandients of the standard basis on T̂, meaning: ∇bₖ(x) = [Pᵢ'(x₂)Pⱼ(x₁),Pᵢ(x₂)Pⱼ'(x₁)] where Pᵢ is the i-th polynomial of Legendre. It returns a vector of functions. 
""" 
function ∇standard_basis(p₁,p₂,p₃)
    n  = compute_dimension(p₁,p₂,p₃)
    ∇B = Vector{Function}(undef,n)
    k = 1
    for i in 0:p₁
        for j in 0:min(p₂,p₃-i)
            ∇B[k] = x->[dnPl(x[1],j,1)*Pl(x[2],i),Pl(x[1],j)*dnPl(x[2],i,1)]
            k   += 1
        end
    end
    return ∇B
end;
∇standard_basis(p₁,p₂) = ∇standard_basis(p₁,p₂,p₂);
∇standard_basis(p₁::T) where T<:Integer    = ∇standard_basis(p₁,p₁);
∇standard_basis(t::Tuple) = ∇standard_basis(t...);


"""

    matrix_F!(F,B,nodes) 

Builds the matrix F, consisting on the evaluations of the standard basis on the boundary nodes. 
"""
function matrix_F(B,nodes)
    nₙ    = size(nodes,2)
    n     = length(B)
    F     = zeros(nₙ,n)
    for k in eachindex(B)
        @inbounds for i in 1:nₙ
            F[i,k] = B[k](nodes[:,i])
        end
    end
    return F
end;

function matrix_F!(F,B,nodes)
    for k in eachindex(B)
        for i in eachindex(eachcol(nodes))
            F[i,k] = B[k](nodes[:,i])
        end
    end
end

""" 

    matrix_C(F) 
    
Computes the matrix that transforms the standard basis into the mixed basis. F is the matrix of evaluations of the standard base at the boundary nodes. The number of degrees of freedom is ``n=n_N+n_B``, where ``n_N`` is the number of nodes at the boundary (and consequently: of nodal functions) and ``n_B`` the number of bubble functions. Hence, the size of F is ``n_N × n``.

"""

function matrix_C(B,nod)
    nₙ = size(nod,2)
    n  = length(B)
    F     = matrix_F(B,nod)
    U,Σ,V = svd!(F,full=true)
    S     = Diagonal(1. ./Σ)
    V₁ = view(V,:,1:nₙ)
    V₂ = view(V,:,nₙ+1:n)
    C  = SMatrix{n,n}([V₁*view(S,:,:)*view(U',:,:) V₂])
end

