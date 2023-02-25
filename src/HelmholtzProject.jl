"""
    HelmholtzProject

Project-package containing  code shared by scripts and notebooks in the project
"""
module HelmholtzProject

using LinearAlgebra
using ExtendableGrids
using ExtendableSparse
using SimplexGridFactory
using Triangulate
using Arpack
using SparseArrays
using KrylovKit
using JacobiDavidson
using GenericArpack

# FEM library, implemented as a  Pluto notebook without package dependencies when included here
include("femlib.jl")

export assemble_reaction_diffusion, Dirichlet, mass_matrix_full, mass_matrix_lumped
export assemble_helmholtz
export fenorms, interpolate

include("eigenlib.jl")
export lapack_eigen, arpack_eigen, krylovkit_eigen, jd_eigen, garpack_eigen

# Example grids,  implemented as a  Pluto notebook without package dependencies when included here
include("grids.jl")
export hminmax, rect_t, rect_r
end # module
