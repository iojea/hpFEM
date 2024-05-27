module hpFEM

using Markdown, LinearAlgebra, StaticArrays, SparseArrays, Dictionaries
using Base.Iterators, Base.Threads
using Triangulate, Printf, ElasticArrays, ExactPredicates, ElasticArrays
using CommonSolve, CuthillMcKee, LegendrePolynomials, GrundmannMoeller
using LinearSolve, Makie, ColorSchemes

import Base: @propagate_inbounds, getindex, isequal, in, hash, show, copy
import StaticArraysCore: check_array_parameters,convert_ntuple


include("hptuple.jl")
include("triangles.jl")
include("edges.jl")
include("meshes.jl")
include("show.jl")
include("refine.jl")
include("plots.jl")
include("basis.jl")
include("matrices.jl")
include("estimates.jl")
include("solvers.jl")

export TupleHP,TriangleHP, EdgeHP, MeshHP, EdgeProperties
export edges, degree, marker, mark!, setdegree!, longestedge
export circmesh, circmesh_graded_center, mark!, estim_distance_origin, rectmesh, l_mesh, l_graded,triangle
export normals,in
export rhs, mass, stiff
export ConstantCoeff, BoundaryConditions, ConstantCoeffProblem
export find_neighbors,normals
export plotmeshhp, plotsolhp

end #module
