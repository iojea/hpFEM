
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
    (;points,trilist,edgelist) = mesh
    # points    = getproperty(mesh,p[:points][]) 
    # trilist = getproperty(mesh,p[:trilist][]) 
    # edgelist     = getproperty(mesh,p[:edgelist][]) 
    noreftris = filter(!ismarked,trilist)
    greentris = filter(isgreen,trilist)
    bluetris  = filter(isblue,trilist)
    redtris   = filter(isred,trilist)
    if !isempty(noreftris)
        noref = hcat([Vector(nodes(t)) for t in triangles(noreftris)]...)'
        poly!(p,points',noref,color=:gray,strokewidth=-0.75,overdraw=false)#,strokecolor=:white,strokewidth=0.25)
    end
    if !isempty(greentris)
        green = hcat([Vector(nodes(t)) for t in triangles(greentris)]...)'
        poly!(p,points',green,color=:olivedrab,strokewidth=-0.75,overdraw=false)#,strokecolor=:lightgray,strokewidth=0.25)
    end
    if !isempty(bluetris)
        blue  = hcat([Vector(nodes(t)) for t in triangles(bluetris)]...)'
        poly!(p,points',blue,color=:dodgerblue,strokewidth=-0.75,overdraw=false)#,strokecolor=:lightgray,strokewidth=0.25)
    end
    if !isempty(redtris)
        red   = hcat([Vector(nodes(t)) for t in triangles(redtris)]...)'
        poly!(p,points',red,color=:tomato,strokewidth=-0.75,overdraw=false)#,strokecolor=:white,strokewidth=0.25)
    end

    for e in edges(edgelist)
        x = Vector(points[1,nodes(e)])
        y = Vector(points[2,nodes(e)])
        lines!(p,x,y,linewidth=p[:linewidth][],color= begin ismarked(edgelist[e]) ? :black : :white end)
    end
    if p[:annotate][]
        for (i,dot) in enumerate(eachcol(points))
            text!(p,Point2f(dot...),text=string(i))
        end
    end
    end
    return p
end


@recipe(PlotSolHP, mesh,u) do scene
    Attributes(
               title = "",
            )
end

function Makie.plot!(p::PlotSolHP)
    (;mesh,u) = p
    lift(p[1]) do mesh
        (;points,trilist,edgelist) = mesh
        lift(p[2]) do u
            minu = minimum(u)
            maxu = maximum(u)
            cr = 0:1;#range(start=minu,stop=maxu,length=3length(u))
            tris = hcat([[t...] for t in triangles(trilist)]...)'
            poly!(p,points',tris,color=(u .-minu)/(maxu-minu),colormap=:coolwarm,colorrange=cr)
        end
    end
    return p
end

function plot_degs(mesh::MeshHP)
    (;points,edgelist) = mesh
    degs = degree.(edgelist)
    nc  = Int(maximum(degs))
    mc  = Int(minimum(degs))
    pal = palette(:blues,nc-mc+1)
    f = Figure()
    Axis(f[1,1])
    for e in edgelist
        nod = nodes(e)
        x = Vector(points[1,nod])
        y = Vector(points[2,nod])
        lines!(x,y,overdraw=true,linewidth=1,color=pal[degree(e)])
    end
    f
end