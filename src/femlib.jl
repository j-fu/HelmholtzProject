### A Pluto.jl notebook ###
# v0.19.6

using Markdown
using InteractiveUtils

# ╔═╡ ea82461b-60a5-4f49-8888-b9fe35cc5a3d
if isdefined(Main, :PlutoRunner)
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    using Revise
    using PlutoUI
    inpluto = true
else
    inpluto = false
end;

# ╔═╡ 72bfdc02-c34c-4a8f-ac65-a158ca33e562
md"""
The contents of this pluto notebook is made available as part of the HelmholtzProject project-package. In that case all pluto specific actions and environment activations are disabled.
"""

# ╔═╡ bde87e19-0ed1-4f24-a453-3693bde5b4d9
if inpluto
    TableOfContents()
end

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
md"""
# P1 finite elements for HelmholtzProject
"""

# ╔═╡ 02d9e61e-271d-4468-acb9-08b7e43f49e3
md"""
## Local stiffness matrix calculation
"""

# ╔═╡ 769773f2-6ebd-4731-a33f-675fc608b33c
md"""
Assume the simplex is given by points  ``a_1 ... a_{d+1}``
Then the barycentric coordinates ``Λ(x)=(λ_1(x) ... λ_{d+1}(x))^T`` are given by
```math
\begin{align}
\sum_{i=1}^{d+1} λ_i(x) &=1\\
\sum_{i=1}^{d+1} a_i λ_i(x)&=x
\end{align}
```
Equivalently,

```math
C Λ = X^+   
```
where ``X^+= (1, x_1 ... x_d)^T`` and ``Λ=(λ_1(x)\dots λ_{d+1}(x))^T``

and ``C`` is the extended matrix of point coordinates
```math
\begin{pmatrix}
  1 & \dots & 1\\
  a_1 & \dots & a_{d+1}
\end{pmatrix}
```

Consequently,
```math
    Λ =C^{-1} X^+   
```

```math 
G=\nabla Λ= \begin{pmatrix}
\ast &\partial_1\lambda_1& \dots& \partial_d \lambda_{1}\\ 
\ast &\vdots &  & \vdots \\
\ast &\partial_1\lambda_{d+1}& \dots& \partial_d \lambda_{d+1} 
\end{pmatrix}
```
-- the matrix of the gradients of  ``λ_j`` can be defined by calculating the partial derivatives of the defining equation:
```math
    G =C^{-1} \mathbb{1}^+   
```
where ``\mathbb{1}^+=\mathrm{diag}(0, 1 \stackrel{n}{...} 1)^T`` is a diagonal matrix and we ignore the first column of the result.

"""

# ╔═╡ 412d7923-4382-4cb9-b228-ff7972bb6af9
"""
     coordmatrix!(C,coord, cellnodes,k)

Assemble extended matrix C of point coordinates.
"""
function coordmatrix!(C, coord, cellnodes, k)
    spacedim = size(coord, 1)
    celldim = size(cellnodes, 1)
    for jj = 1:celldim
        C[1, jj] = 1
        for ii = 1:spacedim
            C[ii + 1, jj] = coord[ii, cellnodes[jj, k]]
        end
    end
end

# ╔═╡ 513991b5-cb2a-439a-8e25-ad84e392d7c3
"""
    gradient!(G,C,factdim,I)

Calculate the gradients gradients 
in rows ``2... d+1`` and  columns ``2 ... d+1`` of G.
"""
function gradient!(G, C, factdim, I)
    clu = lu!(C)
    vol = abs(det(clu)) / factdim
    ldiv!(G, clu, I)
    return vol
end

# ╔═╡ e57e2ef0-8e58-4fe1-89d2-2a9989be2427
"""
    scalpro(G,dim,jl,il)
 Scalar product of columns jl and il of ``d x d`` matrix G, ignoring
 the respective first entries.
"""
function scalpro(G, dim, jl, il)
    s = 0.0
    @inbounds @simd for k = 1:dim
        s += G[jl, k + 1] * G[il, k + 1]
    end
    return s
end

# ╔═╡ 207122e2-11c6-4020-8b19-8d56b9fe9f6b
"""
    stiffness!(S,dim,G,vol)
   Local stiffness matrix for Laplace problem
"""
function stiffness!(S, dim, G)
    @inbounds for il = 1:(dim + 1)
        S[il, il] = scalpro(G, dim, il, il)
        for jl = (il + 1):(dim + 1)
            S[il, jl] = scalpro(G, dim, jl, il)
            S[jl, il] = S[il, jl]
        end
    end
    return S
end

# ╔═╡ f253f9bc-e031-4500-871a-d42b05db9189
"""
    stiffness!(S,dim,permeability_matrix,G)

Local stiffness matrix for anisotropic problem
"""
function stiffness!(S, dim, permeability_matrix, G)
    @inbounds for il = 1:(dim + 1)
        Gil = permeability_matrix * G[il, 2:(dim + 1)]
        S[il, il] = dot(Gil, G[il, 2:(dim + 1)])
        for jl = (il + 1):(dim + 1)
            S[il, jl] = dot(Gil, G[jl, 2:(dim + 1)])
            S[jl, il] = S[il, jl]
        end
    end
    return S
end

# ╔═╡ 7596ac52-029c-48c8-a207-c7f3357e5417
"""
    area(coord, bfacenodes,k, DIM)

Boundary face area
"""

# ╔═╡ b5dcd565-b59c-4c69-a5f4-a1e0e6f0073d
area(coord, bfacenodes, k, ::Type{Val{1}}) = 1

# ╔═╡ e8fac68c-0144-4430-9735-a9736872d967
function area(coord, bfacenodes, k, ::Type{Val{2}})
    dx = coord[1, bfacenodes[2, k]] - coord[1, bfacenodes[1, k]]
    dy = coord[2, bfacenodes[2, k]] - coord[2, bfacenodes[1, k]]
    sqrt(dx^2 + dy^2)
end

# ╔═╡ e0191b9c-d4fd-41ff-860a-259e18691dab
function area(coord, bfacenodes, k, ::Type{Val{3}})
    d21x = coord[1, bfacenodes[2, k]] - coord[1, bfacenodes[1, k]]
    d21y = coord[2, bfacenodes[2, k]] - coord[2, bfacenodes[1, k]]
    d21z = coord[3, bfacenodes[2, k]] - coord[3, bfacenodes[1, k]]

    d31x = coord[1, bfacenodes[3, k]] - coord[1, bfacenodes[1, k]]
    d31y = coord[2, bfacenodes[3, k]] - coord[2, bfacenodes[1, k]]
    d31z = coord[3, bfacenodes[3, k]] - coord[3, bfacenodes[1, k]]

    nx = (d21y * d31z - d31y * d21z)
    ny = (d21z * d31x - d31z * d21x)
    nz = (d21x * d31y - d31x * d21y)

    sqrt(nx^2 + ny^2 + nz^2) / 2
end

# ╔═╡ 7f5be01e-6a05-44ba-9b7e-46b29d5c2ada
md"""
## Local mass matrix calculation
"""

# ╔═╡ 71ee16dd-140c-4709-9e3b-d0e6bd0a3d75
mass_matrix_full(::Type{Val{1}}) = 1 / 6 * [2 1; 1 2.0]

# ╔═╡ d2e6aea7-3035-4d03-b531-44ac98d81c68
mass_matrix_lumped(::Type{Val{1}}) = 1 / 2 * [1 0; 0 1.0]

# ╔═╡ ff0bb725-86b6-479b-a684-4b6d7dc42a4a
mass_matrix_full(::Type{Val{2}}) = 1 / 12 * [2 1 1; 1 2 1; 1 1 2.0]

# ╔═╡ 68fd1efc-0af4-4886-86ca-428ea7ac927e
mass_matrix_lumped(::Type{Val{2}}) = 1 / 3 * [1 0 0; 0 1 0; 0 0 1.0]

# ╔═╡ 468e879d-4dce-46c8-8138-fc74daca6c6d
mass_matrix_full(::Type{Val{3}}) = 5 / 20 * [2 1 1 1; 1 1 2 1; 1 1 2.0 1; 1 1 1 2]

# ╔═╡ 92ebcbdc-642c-417a-9b3f-1e084ad349b1
mass_matrix_lumped(::Type{Val{3}}) = 1 / 4 * [1 0 0 0; 0 1 0 0; 0 0 1.0 0; 0 0 0 1.0]

# ╔═╡ 029d458e-cbe5-45b4-8d9d-faeee1196dd6
md"""
## Assembly for reaction-diffusion
"""

# ╔═╡ 6eb84b8e-bf03-4ac6-84c9-89b190ffe4eb
"""
    Dirichlet()

Marker for Dirichlet boundary conditions
"""
Dirichlet() = 1.0e30

# ╔═╡ 663c5f2d-1e6a-426c-9a74-0c55cd77288f
"""
    assemble_rection_diffusion!(matrix,rhs,grid, k,r,f,b,m)

Assemble matrix and right hand side for P1 finite element discretization of diffusion problem.
"""
function assemble_reaction_diffusion!(A_h, # Global stiffness matrix
                                      F_h, # Right hand side of FEM problem
                                      grid, # Discretization grid 
                                      K::TK, # Diffusion coefficient
                                      R::TR, # Reaction cooefficient
                                      F::TF, # Right hand side
                                      B::TB,  # Boundary function
                                      MMatrix::TM) where {TK, TR, TF, TB, TM}
    coord = grid[Coordinates]
    cellnodes = grid[CellNodes]
    cellregions = grid[CellRegions]
    dim = size(coord, 1)
    nnodes = dim + 1
    factdim::Float64 = factorial(dim)
    S = zeros(nnodes, nnodes) # local stiffness matrix
    C = zeros(nnodes, nnodes)  # local coordinate matrix
    BC = zeros(dim, dim)     # local boundary coordinate matrix
    G = zeros(nnodes, nnodes) # shape function gradients
    I = Matrix(Diagonal(ones(nnodes)))
    M = MMatrix(Val{dim})

    ncells = size(cellnodes, 2)
    for icell = 1:ncells
        iregion = cellregions[icell]
        coordmatrix!(C, coord, cellnodes, icell)
        vol = gradient!(G, C, factdim, I)
        stiffness!(S, dim, G)
        avg_K = 0.0
        avg_R = 0.0
        for il = 1:nnodes
            @views avg_K += K(iregion, coord[:, cellnodes[il, icell]])
            @views avg_R += R(iregion, coord[:, cellnodes[il, icell]])
        end
        avg_K /= (nnodes)
        avg_R /= (nnodes)
        for il = 1:nnodes
            i = cellnodes[il, icell]
            for jl = 1:nnodes
                j = cellnodes[jl, icell]
                A_h[i, j] += vol * (avg_K * S[il, jl] + avg_R * M[il, jl])
            end
            @views F_h[i] += F(iregion, coord[:, cellnodes[il, icell]]) * vol / (nnodes)
        end
    end

    # Boundary part with penalty method
    bfacenodes = grid[BFaceNodes]
    nbfaces = size(bfacenodes, 2)
    bfaceregions = grid[BFaceRegions]
    for ibface = 1:nbfaces
        a = area(coord, bfacenodes, ibface, Val{dim})
        for idim = 1:dim
            i1 = bfacenodes[idim, ibface]
            @views bcfac, bcval = B(bfaceregions[ibface], coord[:, i1])
            if bcfac == Dirichlet()
                A_h[i1, i1] += bcfac
                F_h[i1] += bcfac * bcval
            else # this corresponds to a "lumped" matrix
                A_h[i1, i1] += a * bcfac / dim
                F_h[i1] += a * bcval / dim
            end
        end
    end
end

# ╔═╡ f4205e3c-1a07-493f-aa94-c1e197140be6
md"""
Finite element approximation for 

```math
\begin{aligned}
 -\nabla k \nabla u + ru &= f\quad \text{in}\; \Omega\\
 k\nabla u\cdot \vec n +  \alpha u & = g\quad  \text{on}\; \partial\Omega
\end{aligned}
```
"""

# ╔═╡ 865c31aa-b045-4b8e-a344-c29116dcc194
"""
     A,f=assemble_reaction_diffusion(grid,k,r,f,b)

Assemble reaction-diffusion problem.
Parameters:
- grid: 1/2/3D ExtendableGrid
- k: `(iregion, coord) -> Number` function returning value of diffusion coefficient
- r: `(iregion, coord) -> Number` function returning reaction coefficient
- f: `(iregion, coord) -> Number` function returning value of right hand side
- b: `(ibregion, coord) -> (α,g) function returning coefficient α and boundary value g. If α==Dirichlet(), assume Dirichlet boundary condition `u=g`. 
- m: either `mass_matrix_full` or `mass_matrix_lumped` 
Returns assembled matrix and right hand side
"""
function assemble_reaction_diffusion(grid;
                                     k = (i, coord) -> 1,
                                     r = (i, coord) -> 0,
                                     f = (i, coord) -> 0,
                                     b = (i, coord) -> (Dirichlet(), 0),
                                     m = mass_matrix_full)
    n = num_nodes(grid)
    A_h = ExtendableSparseMatrix(n, n)
    F_h = zeros(n)
    assemble_reaction_diffusion!(A_h, F_h, grid, k, r, f, b, m)
    flush!(A_h)
    A_h.cscmatrix, F_h
end

# ╔═╡ 3fecabf8-456a-4a33-8edd-5c13c8739358
md"""
## Assembly for Helmholtz eigenvalue problem
"""

# ╔═╡ 4dc2030e-207d-4bd2-b718-606d3e0445e4
md"""
```math
-\Delta u - nk^2 u = \lambda u
```
"""

# ╔═╡ fe3d9d51-ac68-486d-b4e6-e49c84d1b616
"""
    assemble_helmholtz!(A,M,grid:k,n,b,m)

Assemble matrix and mass matrix for P1 finite element discretization of
Helmholtz problem
"""
function assemble_helmholtz!(A_h, # Global stiffness matrix
                             M_h, # Global mass matrix
                             grid,
                             k,
                             N::TN,
                             B::TB,
                             MMatrix::TM) where {TN, TB, TM}
    coord = grid[Coordinates]
    cellnodes = grid[CellNodes]
    cellregions = grid[CellRegions]
    dim = size(coord, 1)
    factdim::Float64 = factorial(dim)
    S = zeros(dim + 1, dim + 1)
    C = zeros(dim + 1, dim + 1)
    BC = zeros(dim, dim)
    G = zeros(dim + 1, dim + 1)
    I = Matrix(Diagonal(ones(dim + 1)))
    M = MMatrix(Val{dim})
    ncells = size(cellnodes, 2)
    for icell = 1:ncells
        iregion = cellregions[icell]
        coordmatrix!(C, coord, cellnodes, icell)
        vol = gradient!(G, C, factdim, I)
        stiffness!(S, dim, G)
        avg_N = 0.0
        for il = 1:(dim + 1)
            @views avg_N += N(iregion, coord[:, cellnodes[il, icell]])
        end
        avg_N /= (dim + 1)
        for il = 1:(dim + 1)
            i = cellnodes[il, icell]
            for jl = 1:(dim + 1)
                j = cellnodes[jl, icell]
                A_h[i, j] += vol * (S[il, jl] - k^2 * avg_N * M[il, jl])
                M_h[i, j] += vol * M[il, jl]
            end
        end
    end
    bfacenodes = grid[BFaceNodes]
    nbfaces = size(bfacenodes, 2)
    bfaceregions = grid[BFaceRegions]
    for ibface = 1:nbfaces
        a = area(coord, bfacenodes, ibface, Val{dim})
        for idim = 1:dim
            i1 = bfacenodes[idim, ibface]
            @views bcfac, bcval = B(bfaceregions[ibface], coord[:, i1])
            if bcfac == Dirichlet()
                A_h[i1, i1] += bcfac
            else # this corresponds to a "lumped" matrix
                A_h[i1, i1] += a * bcfac / dim
            end
        end
    end
end

# ╔═╡ ba649254-054d-470c-a61e-7709a00467a3
"""
     A,M=assemble_helmholtz(grid;k,n,b,m)

Assemble Helmholtz problem with homogeneous Neumann boundary conditions
Parameters:
- grid: 1/2/3D ExtendableGrid
- k: complex wave numner
- n: `(iregion, coord) -> Number`: diffraction index
- b: `(iregion, coord) -> Number,Number`: boundary function
- m: either `mass_matrix_full` or `mass_matrix_lumped` 
Returns assembled matrix and mass matrix for generalized eigenvalue problem
"""
function assemble_helmholtz(grid;
                            k = 1.0 + 0im,
                            n = (i, coord) -> 1,
                            b = (i, coord) -> (0, 0),
                            m = mass_matrix_full)
    N = num_nodes(grid)
    A_h = ExtendableSparseMatrix(Complex{Float64}, N, N)
    M_h = ExtendableSparseMatrix(Float64, N, N)
    assemble_helmholtz!(A_h, M_h, grid, k, n, b, m)
    flush!(A_h)
    flush!(M_h)
    A_h.cscmatrix, M_h.cscmatrix
end

# ╔═╡ efd4a240-6acc-4457-859d-be58ade2cb98
function fenorms(u, coord, cellnodes)
    l2norm = 0.0
    h1norm = 0.0
    dim = size(coord, 1)
    nnodes = dim + 1
    factdim::Float64 = factorial(dim)
    S = zeros(nnodes, nnodes) # local stiffness matrix
    C = zeros(nnodes, nnodes)  # local coordinate matrix
    G = zeros(nnodes, nnodes) # shape function gradients
    I = Matrix(Diagonal(ones(nnodes)))
    ncells = size(cellnodes, 2)
    mmatrix = mass_matrix_full(Val{dim})
    for icell = 1:ncells
        coordmatrix!(C, coord, cellnodes, icell)
        vol = gradient!(G, C, factdim, I)
        stiffness!(S, dim, G)
        for i = 1:nnodes
            for j = 1:nnodes
                uij = u[cellnodes[j, icell]] * conj(u[cellnodes[i, icell]])
                l2norm += uij * vol * mmatrix[j, i]
                h1norm += uij * vol * S[j, i]
            end
        end
    end
    return (real(sqrt(l2norm)), real(sqrt(h1norm)))
end

# ╔═╡ 0515aab6-7e50-43ce-adfa-6ab8c453281f
@doc raw"""
	fenorms(u,grid)

Calculate ``L^2`` and ``H^1`` norms of u. 
"""
fenorms(u, grid) = fenorms(u, grid[Coordinates], grid[CellNodes])

# ╔═╡ 1cce4ce0-2afa-4386-8141-8e17c6ba168b
function interpolate!(u_to, grid_to, u_from, grid_from; eps = 1.0e-14)
    coord = grid_to[Coordinates]
    cn_from = grid_from[CellNodes]
    dim = size(coord, 1)
    nnodes_to = size(coord, 2)
    λ = zeros(dim + 1)
    λ0 = zeros(dim + 1)
    icellstart = 1
    cf = CellFinder(grid_from)
    shuffle = dim == 2 ? [3, 1, 2] : [2, 1]
    for inode_to = 1:nnodes_to
        @views icell_from = gFindLocal!(λ, cf, coord[:, inode_to]; icellstart, eps)
        @assert icell_from > 0
        for i = 1:(dim + 1)
            inode_from = cn_from[i, icell_from]
            u_to[inode_to] += λ[shuffle[i]] * u_from[inode_from]
        end
        icell_start = icell_from
    end
    u_to
end

# ╔═╡ 9017567f-2643-4173-b54f-d73d07bc5d93
"""
	interpolate(grid_to, u_from, grid_from;eps=1.0e-14)

Linear interpolation of function `u_from` on grid `grid_from` to `grid_to`.
"""
function interpolate(grid_to, u_from, grid_from; eps = 1.0e-14)
    u_to = zeros(eltype(u_from), num_nodes(grid_to))
    interpolate!(u_to, grid_to, u_from, grid_from; eps)
end

# ╔═╡ Cell order:
# ╟─72bfdc02-c34c-4a8f-ac65-a158ca33e562
# ╠═ea82461b-60a5-4f49-8888-b9fe35cc5a3d
# ╠═bde87e19-0ed1-4f24-a453-3693bde5b4d9
# ╟─60941eaa-1aea-11eb-1277-97b991548781
# ╟─02d9e61e-271d-4468-acb9-08b7e43f49e3
# ╠═769773f2-6ebd-4731-a33f-675fc608b33c
# ╠═412d7923-4382-4cb9-b228-ff7972bb6af9
# ╠═513991b5-cb2a-439a-8e25-ad84e392d7c3
# ╠═e57e2ef0-8e58-4fe1-89d2-2a9989be2427
# ╠═207122e2-11c6-4020-8b19-8d56b9fe9f6b
# ╠═f253f9bc-e031-4500-871a-d42b05db9189
# ╠═7596ac52-029c-48c8-a207-c7f3357e5417
# ╠═b5dcd565-b59c-4c69-a5f4-a1e0e6f0073d
# ╠═e8fac68c-0144-4430-9735-a9736872d967
# ╠═e0191b9c-d4fd-41ff-860a-259e18691dab
# ╟─7f5be01e-6a05-44ba-9b7e-46b29d5c2ada
# ╠═71ee16dd-140c-4709-9e3b-d0e6bd0a3d75
# ╠═d2e6aea7-3035-4d03-b531-44ac98d81c68
# ╠═ff0bb725-86b6-479b-a684-4b6d7dc42a4a
# ╠═68fd1efc-0af4-4886-86ca-428ea7ac927e
# ╠═468e879d-4dce-46c8-8138-fc74daca6c6d
# ╠═92ebcbdc-642c-417a-9b3f-1e084ad349b1
# ╟─029d458e-cbe5-45b4-8d9d-faeee1196dd6
# ╠═6eb84b8e-bf03-4ac6-84c9-89b190ffe4eb
# ╠═663c5f2d-1e6a-426c-9a74-0c55cd77288f
# ╟─f4205e3c-1a07-493f-aa94-c1e197140be6
# ╠═865c31aa-b045-4b8e-a344-c29116dcc194
# ╟─3fecabf8-456a-4a33-8edd-5c13c8739358
# ╟─4dc2030e-207d-4bd2-b718-606d3e0445e4
# ╠═fe3d9d51-ac68-486d-b4e6-e49c84d1b616
# ╠═ba649254-054d-470c-a61e-7709a00467a3
# ╠═efd4a240-6acc-4457-859d-be58ade2cb98
# ╠═0515aab6-7e50-43ce-adfa-6ab8c453281f
# ╠═1cce4ce0-2afa-4386-8141-8e17c6ba168b
# ╠═9017567f-2643-4173-b54f-d73d07bc5d93
