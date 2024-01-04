module FEMhp

using LinearAlgebra
using Dictionaries
using StaticArrays
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

include("triangles.jl")
include("edges.jl")
include("dicthp.jl")
include("meshes.jl")
include("refine.jl")
include("plots.jl")
# include("basis.jl")
# include("matrices.jl")
# include("solvers.jl")

export TriangleHP, EdgeHP, DictHP, MeshHP, circular_mesh, circular_mesh_graded_to_center, mark!, mark_triangles!, estim_distance_origin, sort_degrees, plotmeshhp,set_degree!,degrees_of_freedom, degree
end
