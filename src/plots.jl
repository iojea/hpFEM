
@recipe(PlotMeshHP, mesh) do scene
    Attributes(
               linewidth = 1.5,
               annotate  = false,
               title = "",
            )
end

function Makie.plot!(p::PlotMeshHP)
    (;mesh) = p
    lift(p[1]) do mesh
    (;points,triangles,edges) = mesh
    # points    = getproperty(mesh,p[:points][]) 
    # triangles = getproperty(mesh,p[:triangles][]) 
    # edges     = getproperty(mesh,p[:edges][]) 
    noreftris = filter(!ismarked,triangles)
    greentris = filter(isgreen,triangles)
    bluetris  = filter(isblue,triangles)
    redtris   = filter(isred,triangles)
    if !isempty(noreftris)
        noref = hcat([Vector(t) for t in keys(noreftris)]...)'
        poly!(p,points',noref,color=:gray,strokewidth=-0.75,overdraw=false)#,strokecolor=:white,strokewidth=0.25)
    end
    if !isempty(greentris)
        green = hcat([Vector(t) for t in keys(greentris)]...)'
        poly!(p,points',green,color=:olivedrab,strokewidth=-0.75,overdraw=false)#,strokecolor=:lightgray,strokewidth=0.25)
    end
    if !isempty(bluetris)
        blue  = hcat([Vector(t) for t in keys(bluetris)]...)'
        poly!(p,points',blue,color=:dodgerblue,strokewidth=-0.75,overdraw=false)#,strokecolor=:lightgray,strokewidth=0.25)
    end
    if !isempty(redtris)
        red   = hcat([Vector(t) for t in keys(redtris)]...)'
        poly!(p,points',red,color=:tomato,strokewidth=-0.75,overdraw=false)#,strokecolor=:white,strokewidth=0.25)
    end

    for e in keys(edges)
        x = Vector(points[1,e])
        y = Vector(points[2,e])
        lines!(p,x,y,linewidth=p[:linewidth][],color= begin ismarked(edges[e]) ? :black : :white end)
    end
    if p[:annotate][]
        for (i,dot) in enumerate(eachcol(points))
            text!(p,Point2f(dot...),text=string(i))
        end
    end
    end
    return p
end


function plot_degs(mesh::MeshHP)
    (;points,edges) = mesh
    degs = degree.(edges)
    nc  = Int(maximum(degs))
    mc  = Int(minimum(degs))
    pal = palette(:blues,nc-mc+1)
    f = Figure()
    Axis(f[1,1])
    for e in edges
        nod = nodes(e)
        x = Vector(points[1,nod])
        y = Vector(points[2,nod])
        lines!(x,y,overdraw=true,linewidth=1,color=pal[degree(e)])
    end
    f
end
