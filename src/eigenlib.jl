### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ ea82461b-60a5-4f49-8888-b9fe35cc5a3d
if isdefined(Main,:PlutoRunner)
    using Pkg
	Pkg.activate(joinpath(@__DIR__,"..")) 
    using Revise
	using PlutoUI
	using Arpack
	using GenericArpack
	using ExtendableSparse
	using SparseArrays
	using LinearAlgebra
	using KrylovKit
	using JacobiDavidson
	inpluto=true
else
	inpluto=false
end;

# ╔═╡ 72bfdc02-c34c-4a8f-ac65-a158ca33e562
md"""
The contents of this pluto notebook is made available as part of the HelmholtzProject project-package. In that case all pluto specific actions and environment activations are disabled.
"""

# ╔═╡ bde87e19-0ed1-4f24-a453-3693bde5b4d9
if inpluto TableOfContents() end  	

# ╔═╡ 27699fe1-486c-43ac-aa96-07fdef4cb9be
function normalize!(v,M)
	nev=size(v,2)
	for i=1:nev
		@views nm=sqrt(dot(M*v[:,i],v[:,i]))
		if abs(nm)>0 
		@views v[:,i]./=nm 
		end
	end
	v
end

# ╔═╡ 7fc7eaf4-dbcc-4f67-ab01-c46a201a3e13
function lapack_eigen(A,M;nev=10,maxiter=300,tol=0.0)
		e=eigen(Matrix(A),Matrix(M))
	e.values, normalize!(e.vectors,M)
end

# ╔═╡ c625ff6e-1d6f-4620-bcb0-45b3ca1b1762
function arpack_eigen(A,M,;nev=10,maxiter=300,tol=1.0e-10)
	λ,v=Arpack.eigs(Symmetric(A),M,ritzvec=true,which=:SR;nev,maxiter,tol)
	λ,normalize!(v,M)
end

# ╔═╡ 21992eed-b614-41b2-9bd8-54d3d0e5c5fd
function garpack_eigen(A,M,;nev=10,maxiter=300,tol=1.0e-10)
	λ,v=GenericArpack.symmeigs(Symmetric(A),Symmetric(M),nev;maxiter,which=:SM,tol)
	λ,normalize!(v,M)
end

# ╔═╡ 064ece9e-23d3-4ad6-94a1-250df864e071
function krylovkit_eigen(A,M; nev=10,maxiter=300, tol=1.0e-10)
λ, v, info = geneigsolve((Symmetric(A),Symmetric(M));
	which=:SR,
	howmany=nev,
    tol)
λ,normalize!(v,M)
end

# ╔═╡ cee2bc6b-53f8-482c-af53-9ff246fcc003
function jd_eigen(A,M; nev=10,maxiter=300, tol=1.0e-9)
	global verbose=false
	target=Near(-100_000+0.1im)

	pschur, residuals = jdqz(A, M,	
    pairs = nev,
    target= target,
    tolerance = tol,
    subspace_dimensions = 10:nev+10,
    max_iter = maxiter,
    verbosity=1
    )
    v=copy(pschur.Q.all)
	pschur.alphas./pschur.betas, normalize!(v,M)
end

# ╔═╡ Cell order:
# ╟─72bfdc02-c34c-4a8f-ac65-a158ca33e562
# ╠═ea82461b-60a5-4f49-8888-b9fe35cc5a3d
# ╠═bde87e19-0ed1-4f24-a453-3693bde5b4d9
# ╠═27699fe1-486c-43ac-aa96-07fdef4cb9be
# ╠═7fc7eaf4-dbcc-4f67-ab01-c46a201a3e13
# ╠═c625ff6e-1d6f-4620-bcb0-45b3ca1b1762
# ╠═21992eed-b614-41b2-9bd8-54d3d0e5c5fd
# ╠═064ece9e-23d3-4ad6-94a1-250df864e071
# ╠═cee2bc6b-53f8-482c-af53-9ff246fcc003
