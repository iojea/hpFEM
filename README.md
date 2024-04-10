
<p align="center">
    <img width="200px" src="https://github.com/iojea/hpFEM/blob/main/docs/assets/logo_text.png"/>
</p> 

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://iojea.github.io/FEMhp.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://iojea.github.io/FEMhp.jl/dev/)
<!--[![Build Status](https://github.com/iojea/FEMhp.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/iojea/FEMhp.jl/actions/workflows/CI.yml?query=branch%3Amain)-->
[![Build Status](https://travis-ci.com/iojea/FEMhp.jl.svg?branch=main)](https://travis-ci.com/iojea/FEMhp.jl)

A pure Julia library for solving 2D partial differential equations using a posteriori $hp$-refinement on triangular meshes. 

***hpFEM*** is under heavy development. At the moment it implements its own $h$ refinement algorithm, following a *red-blue-green* marking strategy and a solver for problems of the form 
$$-\alpha\Delta u + \beta u = f,$$
with homogeneous Dirichlet boundary conditions. An algorithm for $p$ refinement is also implemented, but it is not linked to the refiner, yet.  

The main idea is as follows: each edge $e$ of the triangulation is assigned a degree $p_e$. On each triangle $T$ the solution is approximated using a polynomial $P$ of degree $\max_{e\subset T} p_e$ such that $\textrm{deg}(P|_e) = p_e$. If the degrees assigned to the edges of a triangle are $p_1\le p_2\le p_3$ the restriction $p_3\le p_1 + p_2$ is necessary for the polynomial space to be well-defined.
 
Currenly, it is possible to generate meshes for arbitrary polygons and manually assign degrees to each edge. A $p$-conformity routine is run in order to make the assignment feasible. It is also possible to perform $h$-refinement of the mesh. Triangulate(https://github.com/JuliaGeometry/Triangulate.jl), a Julia wrapper for the Triangle(https://www.cs.cmu.edu/~quake/triangle.html) mesh generator is used for the creation of the initial mesh. 
