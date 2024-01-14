module FEMhp

using Markdown
using LinearAlgebra
using Dictionaries
using StaticArrays
using SparseArrays
using Base.Iterators
using Makie
using ColorSchemes
using Triangulate
using Printf
using ElasticArrays 
using LegendrePolynomials
using GrundmannMoeller
using Base.Threads
using ExactPredicates
using LinearSolve
using CuthillMcKee

include("hptuple.jl")
include("triangles.jl")
include("edges.jl")
include("meshes.jl")
include("refine.jl")
include("plots.jl")
 include("basis.jl")
 include("matrices.jl")
# include("solvers.jl")

export TriangleHP, EdgeHP, MeshHP
export edges, nodes, degree, marker, mark!, setdegree!, longestedge
export circmesh, circular_mesh_graded_to_center, mark_triangles!, estim_distance_origin, rectmesh
export compute_rhs, compute_mass, compute_stiff_std
export plotmeshhp, degrees_of_freedom, degrees_of_freedom_by_edge
end
