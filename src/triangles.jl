"""
    TriangleHP{I<:Integer,P<:Integer,Bool}

A triangle of an HPMesh, formed by: 

+ `vertices::Tuple{I,I,I}`: indicating the indices of the vertices that form the triangle.
+ `props::TriangleProps{P,Bool}`: with properties of the triangle. 
"""
struct TriangleHP{I<:Integer,P<:Integer}
    nodes::SVector{3,I}
    refine::Ref{P}
    #props::TriangleProps{P,Bool}
end


"""
    TriangleHP(v)

Builds a triangle with vertices defined in the `Tuple` or `Array` `v`, longest edge defined in the `Tuple` `l`, refinement indication in the 

"""
TriangleHP(v::V,r::P) where {V<:AbstractArray,P<:Integer} = TriangleHP(SVector{3}(v),Ref(r))
TriangleHP(v::V) where V<:AbstractArray = TriangleHP(v,0)
function triangle(v::V,p::Matrix{F}) where {V<:Union{Tuple,AbstractArray},F<:AbstractFloat}
    maxi = argmax(norm(p[:,v[mod1(i+1,3)]]-p[:,v[i]]) for i in 1:3)
    TriangleHP(v[mod1.(maxi:maxi+2,3)],UInt8(0))
end


@inline _eval(t::TriangleHP,k)       = t.nodes[mod1(k,3)]
@inline _eval(t::TriangleHP,ind::AbstractRange) = [_eval(t,i) for i in ind]
@inline _eval(t::TriangleHP,v::Array) = [_eval(t,i) for i in v]

#@inline _pfirst(t::TriangleHP)       = t.pfirst
@inline nodes(t::TriangleHP)       = t.nodes
@inline get_edges(t::TriangleHP)   = [SVector{2}(_eval(t,i),_eval(t,i+1)) for i in 1:3]


@inline ismarked(t::TriangleHP)    = t.refine[] > 0
@inline isgreen(t::TriangleHP)     = t.refine[] == 1
@inline isblue(t::TriangleHP)      = t.refine[] == 2
@inline isred(t::TriangleHP)       = t.refine[] == 3
@inline mark!(t::TriangleHP,k)     = t.refine[] = k
@inline mark!(t::TriangleHP)       = mark!(t,3)
@inline longestedge(t::TriangleHP) = SVector{2}(_eval(t,1:2))

#get_psorted_nodes(t::TriangleHP) = _isporiented(t) ? _eval(t,_pfirst(t):_pfirst(t)+2) : _eval(t,_pfirst(t):-1:_pfirst(t)-2)



@inline Base.hash(t::TriangleHP) = Base.hash(t.nodes)
#@inline Base.isequal(t1::TriangleHP,t2::TriangleHP) = Base.isequal(t1.nodes,t2.nodes) 
#@inline Base.isequal(t::SVector{3,I},tri::TriangleHP) where I<:Integer = Base.isequal(t,tri.nodes)
#@inline Base.isequal(tri::TriangleHP,t::SVector{3,I}) where I<:Integer = Base.isequal(t,tri)

function Base.show(io::IO,t::TriangleHP)
    compact = get(io,:compact,false)
    if !compact
        print(io,"(",t.nodes[1],",",t.nodes[2],",",t.nodes[3],")")
        if t.refine == 1
            print(io," => refine green")
        elseif t.refine==2
            print(io," => refine blue")
        elseif t.refine==3
            print(io,"=> refine red")
        end
    else
        show(io,t.nodes)
    end
end

