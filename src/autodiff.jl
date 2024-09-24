struct FEScalarField{M<:MeshHP,F<:AbstractFloat}
    mesh::M
    vals::F
end



function ∇(u,x)
    g = similar(x)
    for i in eachindex(x)
        hx    = Hyper.(x,0.,0.,0.)
        hx[i] = Hyper(x[i],1.,1.,0.)
        g[i]  = ε₁part(u(hx))
    end
    return g
end
∇(u) = x->∇(u,x)

function div(u,x)
    d = 0.
    for i in eachindex(x)
        hx = Hyper.(Float64.(x),0.,0.,0.)
        hx[i] = Hyper(Float64(x[i]),1.,1.,0.)
        d += ε₁part((t->u(t)[i])(hx))
    end
    d
end
div(u) = x->div(u,x)


function Δ(u,x)
    d = 0.
    for i in eachindex(x)
        hx    = Hyper.(Float64.(x),0.,0.,0.)
        hx[i] = Hyper(Float64(x[i]),1.,1.,0.)
        d    += ε₁ε₂part((t->u(t)[i])(hx))
    end
    d
end
Δu(u) = x->Δ(u,x)

