# const PairDict{T,U} = Dictionaries.PairDictionary{T,U,Dictionary{T,U}}
# const EdgeList{I,P} = PairDict{EdgeHP{I},EdgeProperties{P,Bool}} where {I<:Integer,P<:Integer}
# const TriangleList{I,P} = PairDict{TriangleHP{I},TriangleProperties{P}} where {I<:Integer,P<:Integer}

const TriangleList{I,P} = Dictionary{TriangleHP{I},TriangleProperties{P}} where {I,P}
const EdgeList{I,P} = Dictionary{EdgeHP{I},EdgeProperties{P,Bool}} where {I,P}


"""
    MeshHP{F<:AbstractFloat,I<:Integer,P<:Integer}()

A mesh for `HP` finite element methods. Its fields are: 
    + `points::Matrix{F}`: a matrix of points
    + `trilist::FESet{TriangleHP{I,P,B}}`: set of trilist
    + `edgelist::FESet{EdgeHP{I,P,B}}`: set of edgelist

    MeshHP(tri::TriangulationIO)
builds an `MeshHP` from a `Triangulate.TriangulatioIO` struct.  
"""
struct MeshHP{F<:AbstractFloat,I<:Integer,P<:Integer}
    points::ElasticMatrix{F,Vector{F}}
    trilist::TriangleList{I,P}
    edgelist::EdgeList{I,P}    
end
MeshHP(mat::Matrix,tris,edgs) = MeshHP(ElasticMatrix(mat),tris,edgs)
Base.copy(mesh::MeshHP) = MeshHP(deepcopy(mesh.points),deepcopy(mesh.trilist),deepcopy(mesh.edgelist))

function MeshHP(tri::TriangulateIO)
    (;pointlist,trianglelist,edgelist,edgemarkerlist) = tri
    I = eltype(trianglelist)
    P = Int8
    trilist = TriangleList{I,P}([triangle(t,pointlist) for t in eachcol(trianglelist)],[TriangleProperties{P}(0) for _ in 1:size(trianglelist,2)])
    edgelist = dictionary([EdgeHP{I}(e) => EdgeProperties{P,Bool}(1,edgemarkerlist[i],false) for (i,e) in enumerate(eachcol(edgelist))] )
    MeshHP(pointlist,trilist,edgelist)
end

@inline edges(list::T) where T<:EdgeList = keys(list)
@inline triangles(list::T) where T<:TriangleList = keys(list)
@inline edges(mesh::T) where T<:MeshHP = keys(mesh.edgelist)
@inline triangles(mesh::T) where T<:MeshHP = keys(mesh.trilist) 


function same_order(e::EdgeHP{I},elist::EdgeList{I}) where I
    _,(_,k) = gettoken(elist,e)
    oe      = gettokenvalue(keys(elist),k)
    nodes(oe) == nodes(e)
end

function Base.show(io::IO,mime::MIME"text/plain",mesh::MeshHP{I,P}) where {I,P}
    println(io,typeof(mesh))
    header = Markdown.parse("""
        + $(size(mesh.points,2)) nodes.
        + $(length(mesh.trilist)) triangles.
        + $(length(mesh.edgelist)) edges.
    """)
    show(io,mime,header)
    show(io,mime,mesh.points)
    println(io)
    show(io,mime,mesh.trilist)
    println(io)
    show(io,mime,mesh.edgelist)
end
function Base.show(io::IO,mesh::MeshHP)
    println(io,)
    show(io,mesh.points)
    println(io)
    show(io,mesh.trilist)
    println(io)
    show(io,mesh.edgelist)
end

function psortperm(v)
    i    = argmin(v)
    if v[i] ≤ v[mod1(i+1,3)] ≤ v[mod1(i+2,3)]
        ind = mod1.(i:i+2,3)
    else
        ind = mod1.(i:-1:i-2,3)
    end
    return ind
end

function pedges(t::TriangleHP{I},mesh::MeshHP{F,I,P}) where {F,I,P}
    (;edgelist) = mesh
    eds  = edges(t)
    p    = degree.(getindices(edgelist,eds))
    pind = psortperm(p)
    DegTuple(p[pind]),eds[pind] 
end

function pnodes(t::TriangleHP{I},mesh::MeshHP{F,I,P}) where {F,I,P}
    (;edgelist) = mesh
    eds = edges(t)
    p   = degree.(getindices(edgelist,eds))
    pind = psortperm(p)
    DegTuple(p[pind]),first.(eds[pind])
end

### SOME MESHES


function circmesh(c,r,h)
    n      = Int(1+2π*r÷h)
    θ      = range(start=0,stop=2π,length=n)[1:end-1]
    points = hcat([[c[1]+cos(t),c[2]+sin(t)] for t in θ]...)
    edgelist  = hcat([[i,i+1] for i in 1:length(points)-1]...,
                          [size(points,2),1])
    marks  = ones(length(edgelist))
    tri = TriangulateIO(;pointlist=points,edgelist=edgelist,edgemarkerlist=marks)
    maxa = Printf.@sprintf "%.15f" h^2/2
    angle= Printf.@sprintf "%.15f" 30.
    (tri,) = triangulate("ea$(maxa)q$(angle)Q",tri)
    return MeshHP(tri)
end;
circmesh(r,h) = circmesh(zeros(2),r,h)
circmesh(h)   = circmesh(1,h)


function correct_boundary_circular(mesh::MeshHP)
    (;points,edgelist) = mesh
    for e in keys(edgelist)
        if marker(edgelist[e])==1
            i,j = nodes(e)
            points[:,i] .= points[:,i]/norm(points[:,i])
            points[:,j] .= points[:,j]/norm(points[:,j])
        end
    end
end

function circmesh_graded_center(h,μ;maxiter=4,rec=false)
    k = 0
    mesh = circmesh(h)
    flag = true
    mshs = [copy(mesh)]
    while k≤maxiter && flag
        mark!(estim_distance_origin,mesh,h=h,μ=μ)
        rec ? push!(mshs,copy(mesh)) : nothing
        if count(ismarked,mesh.trilist) > 0
            refine!(mesh)
            correct_boundary_circular(mesh)
            k += 1
            rec ? push!(mshs,copy(mesh)) : nothing
        else
            flag = false
        end
    end
    return rec ? mshs : mesh
end

function rectmesh(a,b,c,d,h)
    x = range(start=a,stop=b,length=Int(1+(b-a)÷h))
    n = length(x)
    y = range(start=c,stop=d,length=Int(1+(d-c)÷h))   
    m = length(y) 
    points   =  [x fill(c,n);
                 fill(b,m-2) y[2:end-1];
                 reverse(x) fill(d,n);
                 fill(a,m-2) reverse(y[2:end-1])]'
    edgelist = hcat([[i,i+1] for i in 1:size(points,2)-1]...,
                          [size(points,2),1])
    marks    = ones(Int8,length(edgelist))
    tri = TriangulateIO(;pointlist=points,edgelist=edgelist,edgemarkerlist=marks)
    maxa = Printf.@sprintf "%.15f" h^2/2
    angle= Printf.@sprintf "%.15f" 30.
    (tri,_) = triangulate("ea$(maxa)q$(angle)Q",tri)
    MeshHP(tri)
end
rectmesh(b,d,h) = rectmesh(0,b,0,d,h)
rectmesh(b,h)   = rectmesh(0,b,0,b,h)


function boundarynodes(mesh::MeshHP{F,I,P}) where {F,I,P}
    (;edgelist) = mesh
    v = zeros(I,sum(==(1),marker.(edgelist)))
    i = 1
    for e in edges(edgelist)
        if marker(edgelist[e])==1
            for j in e
                if j ∉ v
                    v[i] = j
                    i   += 1
                end
            end
        end
    end
    v
end
