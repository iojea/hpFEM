const TriangleHP{I} = TupleHP{3,I}

"""
    TriangleHP{I}(x::Tuple)
    TriangleHP{I}(x1,x2)

constructs an edge formed by a pair of indices of type `I<:Integer`. The type can be inferred from the data:
    
    TriangleHP(x::Tuple)
    TriangleHP(x1,x2,x3)

`TriangleHP{I}` is just an alias for `TupleHP{3,I}`.

"""
TriangleHP(x) = TriangleHP{eltype(x)}(x)


"""
    _eval(t::TriangleHP,k)

Returns the index stored in the triangle at index `k` mod 3. 
""" 
@inline _eval(t::TriangleHP,k) = t[mod1(k,3)]


"""
    edges(t::TriangleHP)

Return a vector of edges with type `EdgeHP`, containing the edges of `t`.
"""
@inline edges(t::TriangleHP)   = [EdgeHP(_eval(t,i),_eval(t,i+1)) for i in 1:3]


""" 
    longestedge(t::TriangleHP)

Returns an `EdgeHP` with the longest edge of `t`.
"""
 @inline longestedge(t::TriangleHP) = EdgeHP(_eval(t,1),_eval(t,2))


"""
  triangle(t::T,p::AbstractMatrix)

constructs a `TriangleHP` where the vertices are given by the columns of `p[:,t]`. Hence, the triangle is defined by the indices stored in `t`, but sorted in such a way that the first edge is the longest. 
"""
function triangle(t::T,p::AbstractMatrix) where {T<:Union{Tuple,AbstractArray}}
    maxi = argmax(norm(p[:,t[mod1(i+1,3)]]-p[:,t[i]]) for i in 1:3)
    TriangleHP{eltype(t)}(t[mod1.(maxi:maxi+2,3)])
end


# Properties

"""
    TriangleProperties{P,F}(refine,η,ηₚ) where {P<:Integer,F<:AbstractFloat}

constructs a `struct` for storing attributes of a triangle. These attributes are:
+ `refine`: 
    - `0`: not marked for refinement. 
    - `1`: marked for refinement of _green_ type.
    - `2`: marked for refinement of _blue_ type.
    - `3`: marked for refinement of _red_ type. 
+ `η`: estimate for the local error. 
+ `ηₚ`: predictor of local error, based on previos estimations.  
The types can be inferred from the data:

    TriangleProperties(refine,η,ηₚ)

If only the `refine` argument is passed, `η` and `ηₚ` are initialized as `0.`
If no arguments are passed, `refine` is initialized as `Int8(0)`.
"""
struct TriangleProperties{P<:Integer,F<:AbstractFloat} 
    refine::Ref{P}
    η::Ref{F}
    ηₚ::Ref{F}
end
TriangleProperties{P,F}(val,η,ηₚ) where {P,F} = TriangleProperties(Ref(P(val)),Ref(F(η)),Ref(F(ηₚ)))
TriangleProperties{P}(val,η,ηₚ) where {P} = TriangleProperties(Ref(P(val)),Ref(η),Ref(ηₚ))
TriangleProperties(val::P,η::F,ηₚ::F) where {P,F} = TriangleProperties(Ref(val),Ref(η),Ref(ηₚ))
TriangleProperties{P}() where {P} = TriangleProperties(P(0),0.,0.)
TriangleProperties() = TriangleProperties(Int8(0),0.,0.)
TriangleProperties(tp::TriangleProperties) = TriangleProperties(tp.refine,tp.η,tp.ηₚ)
TriangleProperties{P}(tp::TriangleProperties) where {P} = TriangleProperties(P(tp.refine[]),tp.η,tp.ηₚ)
TriangleProperties{P,F}(tp::TriangleProperties) where {P,F}= TriangleProperties(P(tp.refine[]),F(tp.η),F(tp.ηₚ))

@inline ismarked(t::TriangleProperties)    = t.refine[] > 0
@inline isgreen(t::TriangleProperties)     = t.refine[] == 1
@inline isblue(t::TriangleProperties)      = t.refine[] == 2
@inline isred(t::TriangleProperties)       = t.refine[] == 3
@inline mark!(t::TriangleProperties,k)     = t.refine[] = k
@inline mark!(t::TriangleProperties)       = mark!(t,3)
@inline setη!(t::TriangleProperties,η)     = t.η[] = η
@inline setηₚ!(t::TriangleProperties,ηₚ)   = t.ηₚ[] = ηₚ


