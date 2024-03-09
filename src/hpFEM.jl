module hpFEM

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
using CommonSolve
include("hptuple.jl")
include("triangles.jl")
include("edges.jl")
include("meshes.jl")
include("show.jl")
include("refine.jl")
include("plots.jl")
include("basis.jl")
include("matrices.jl")
include("solvers.jl")

export TriangleHP, EdgeHP, MeshHP
export edges, nodes, degree, marker, mark!, setdegree!, longestedge
export circmesh, circmesh_graded_center, mark!, estim_distance_origin, rectmesh
export rhs, mass, stiff
export ConstantCoeff, BoundaryConditions, ConstantCoeffProblem
export plotmeshhp, plotsolhp

end
