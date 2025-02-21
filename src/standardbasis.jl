struct StandardBasis{T<:Number,I<:Integer,N}
    x::SVector{2,T}
    p₁::I
    p₂::I
    p₃::I
    lix::LegendreVectorIterator{N,T}
    liy::LegendreVectorIterator{N,T}
    function StandardBasis{T,I}(x,p₁,p₂,p₃) where {T,I}
        -1 .<= x .<= 1 || throw(DomainError("x must lie between -1 and 1."))
        lix = LegendreVectorIterator(view(x,1,:))
        liy = LegendreVectorIterator(view(x,2,:))
        new{T,I}(SVector{2,T}(convert.(T,x)),p₁,p₂,p₃,lix,liy)
    end
end





@inline function iterate(basis::StandardBasis{T,I,N}) where {T,I,N}
    px,stx = iterate(basis.lix)
    py,sty = iterate(basis.liy)    
    st     = (zero(T),zero(T))
    p      = px*py
    return p,(st,stx,sty)
end


function iterate(basis::StandardBasis,state)
            
end
    for i in 0:p₁
        for j in 0:min(p₂,p₃-i)

