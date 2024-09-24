# Show TriangleHP
# Base.show(io::IO,mime::MIME"text/plain",t::TriangleHP) = show(io,mime,Int.(t))
# Base.show(io::IO,t::TriangleHP) = show(io,Int.(t))

#Show TriangleProperties
function Base.show(io::IO, mime::MIME"text/plain", t::TriangleProperties)
    symbs = Dictionary(0:3, [:noref, :green, :blue, :red])
    show(io, symbs[t.refine[]])
end
function Base.show(io::IO, t::TriangleProperties)
    symbs = Dictionary(0:3, [:noref, :green, :blue, :red])
    show(io, symbs[t.refine[]])
end

# Show EdgeHP
# Base.show(io::IO,mime::MIME"text/plain",t::EdgeHP) = show(io,mime,Int.(t))
# Base.show(io::IO,t::EdgeHP) = show(io,Int.(t.nodes[:]))

# Show EdgeProperties
function Base.show(io::IO, mime::MIME"text/plain", t::EdgeProperties)
    if t.marker == 0
        m = :Ω°
    elseif t.marker == 1
        m = :∂𝔇
    elseif t.marker == 2
        m = :∂𝔑
    else
        m = t.marker
    end
    show(io, mime, (t.degree[], m, t.refine[]))
end
function Base.show(io::IO, t::EdgeProperties)
    if t.marker == 0
        m = :Ω°
    elseif t.marker == 1
        m = :Γ𝔇
    elseif t.marker == 2
        m = :Γ𝔫
    else
        m = t.marker
    end
    if t.refine[]
        r = :refine
    else
        r = :noref
    end
    show(io, (t.degree[], m, r))
end

# Show HPList short.
function Base.show(io::IO, d::Union{TriangleList,EdgeList})
    limit = get(io, :limit, false) ? Int64(10) : typemax(Int64)
    comma = false
    print(io, "{")
    for i in keys(d)
        if comma
            print(io, ", ")
        end
        if limit == 0
            print(io, "…")
            break
        end
        show(io, i)
        print(io, " = ")
        show(io, d[i])
        comma = true
        limit -= 1
    end
    print(io, "}")
end


# The "display"-like rendering
function Base.show(io::IO, ::MIME"text/plain", d::Union{TriangleList,EdgeList})
    if isempty(d)
        print(io, "0-element $(typeof(d))")
        return
    end

    # Designed to be efficient for very large sets of unknown lengths

    n_lines = max(
        Int64(3),
        get(io, :limit, false) ? Int64(displaysize(io)[1] - 4) : typemax(Int64),
    )
    n_cols = max(
        Int64(8),
        get(io, :limit, false) ? Int64(displaysize(io)[2] - 4) : typemax(Int64),
    )
    n_lines_top = n_lines ÷ Int64(2)
    n_lines_bottom = n_lines - n_lines_top

    # First we collect strings of all the relevant elements
    top_ind_strs = Vector{String}()
    top_val_strs = Vector{String}()
    bottom_val_strs = Vector{String}()
    bottom_ind_strs = Vector{String}()

    top_lines = Int64(1)
    top_full = false
    top_last_index = Base.RefValue{keytype(d)}()
    for i in keys(d)
        push!(top_ind_strs, sprint(show, i, context = io, sizehint = 0))
        push!(top_val_strs, sprint(show, d[i], context = io, sizehint = 0))
        top_lines += 1
        if top_lines > n_lines_top
            top_full = true
            top_last_index[] = i
            break
        end
    end

    bottom_lines = Int64(1)
    bottom_full = false
    if top_full
        for i in Iterators.reverse(keys(d))
            if bottom_full
                if isequal(i, top_last_index[])
                    bottom_full = false # overide this, we don't need the ⋮
                else
                    bottom_ind_strs[end] = "⋮"
                    bottom_val_strs[end] = "⋮"
                end
                break
            end

            if isequal(i, top_last_index[])
                # Already rendered, we are finished
                break
            end

            push!(bottom_ind_strs, sprint(show, i, context = io, sizehint = 0))
            push!(bottom_val_strs, sprint(show, d[i], context = io, sizehint = 0))

            bottom_lines += 1
            if bottom_lines > n_lines_bottom
                bottom_full = true
                # We check the next element to see if this one should be a ⋮
            end
        end
        ind_strs = vcat(top_ind_strs, reverse(bottom_ind_strs))
        val_strs = vcat(top_val_strs, reverse(bottom_val_strs))
    else
        ind_strs = top_ind_strs
        val_strs = top_val_strs
    end

    if Base.IteratorSize(d) === Base.SizeUnknown()
        if bottom_full
            print(io, "Greater than $(length(ind_strs))-element $(typeof(d))")
        else
            print(io, "$(length(ind_strs))-element $(typeof(d))")
        end
    else
        print(io, "$(length(d))-element $(typeof(d))")
    end

    # Now find padding sizes
    max_ind_width = maximum(textwidth, ind_strs)
    max_val_width = maximum(textwidth, val_strs)
    if max_ind_width + max_val_width > n_cols
        val_width = max_val_width
        ind_width = max_ind_width
        while ind_width + val_width > n_cols
            if ind_width > val_width
                ind_width -= 1
            else
                val_width -= 1
            end
        end
        if ind_width != max_ind_width
            shrink_to!(ind_strs, ind_width)
        end
        if val_width != max_val_width
            shrink_to!(val_strs, val_width)
        end
    else
        ind_width = max_ind_width
    end

    for (ind_str, val_str) in zip(ind_strs, val_strs)
        print(io, "\n ")
        print(io, " "^max(0, ind_width - textwidth(ind_str)))
        print(io, ind_str)
        print(io, " │ ")
        print(io, val_str)
    end
end

function shrink_to!(strs, width)
    for i in keys(strs)
        str = strs[i]
        if textwidth(str) > width
            new_str = ""
            w = 0
            for c in str
                new_w = textwidth(c)
                if new_w + w < width
                    new_str = new_str * c
                    w += new_w
                else
                    new_str = new_str * "…"
                    break
                end
            end
            strs[i] = new_str
        end
    end
end

# Show MeshHP
function Base.show(io::IO, mime::MIME"text/plain", mesh::MeshHP{I,P}) where {I,P}
    println(io, typeof(mesh))
    header = Markdown.parse("""
        + $(size(mesh.points,2)) nodes.
        + $(length(mesh.trilist)) triangles.
        + $(length(mesh.edgelist)) edges.
    """)
    show(io, mime, header)
    println(io)
    show(io, mime, mesh.points)
    println(io)
    println(io)
    show(io, mime, mesh.trilist)
    println(io)
    println(io)
    println(io)
    show(io, mime, mesh.edgelist)
end


function Base.show(io::IO, mesh::MeshHP)
    println(io)
    show(io, mesh.points)
    println(io)
    show(io, mesh.trilist)
    println(io)
    show(io, mesh.edgelist)
end
