## Meshes

The mesh struct is `MeshHP`. It is formed by an `ElasticMatrix` with the coordinates of the vertices, a set of triangles and a set of edges. This sets are given by structs of type `FESet`, which is just a wrapper around an `Indices` struct (from `Dictionaries.jl`). The goal of this wrapper is to be able to define special methods for some functions without interfere with the normal behaviour of `Indices`. 

### Triangles

The struct `TriangleHP` has a field `nodes::SVector{3,I} where I<:Integer` and a field of `TriangleProps`, that contains the attributes of the triangle. `nodes` is always positively oriented. It is a `StaticVector` and not a `Tuple` because it used for indexing over the matrix of vertices. `TriangleProps` is meant to be a hidden struct. It contains: the index of first node of its longest edge, a number indicating if the triangle is marked for refinement `TriangleProps` and a pair of fields, `pfirst` and `poriented` that indicates how the triangle should be read so that the corresponding edge have non-decreasing associated polynomial degrees.  Special methods are defined in order to retrieve this infomation directly from the triangle. 

### Edges
`EdgeHP` is similar to `TriangleHP`
