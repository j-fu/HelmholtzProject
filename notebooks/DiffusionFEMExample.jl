### A Pluto.jl notebook ###
# v0.18.4

using Markdown
using InteractiveUtils

# ╔═╡ f9a43662-326a-46b4-a638-39498fb59ce6
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise
	using ExtendableGrids
	using PlutoVista
	using PlutoUI
	using GridVisualize
	using ExtendableSparse
	using HelmholtzProject
	default_plotter!(PlutoVista)	
end

# ╔═╡ 0f649cee-b670-4c66-9c6f-72e6a93f52a0
TableOfContents(title="")

# ╔═╡ 9eaaa834-0b7b-46a6-a300-a70a5c7af9db
md"""
### 1D Example
"""

# ╔═╡ f514b55e-92b6-418b-8a29-cf0c425453c8
X1=0:0.1:1

# ╔═╡ b1e6e0d8-c22d-453d-9deb-cbde27656689
g1=simplexgrid(X1)

# ╔═╡ 1d888171-2b58-4370-8bf2-63905a3d37f9
gridplot(g1,legend=:lt,resolution=(500,150))

# ╔═╡ 4de4e9c5-ee40-40be-9bed-9f538d45674e
A1,f1=assemble_reaction_diffusion(g1,f=(i,x)->2)

# ╔═╡ 9888e062-57f7-434a-9176-0a3cd8230448
sol1=A1\f1

# ╔═╡ 37f9ca73-d1e6-4afc-b2a3-858a80faa2e5
scalarplot(g1,sol1,resolution=(500,150),label="f",legend=:rt,xlabel="x",ylabel="y")

# ╔═╡ f65baba2-6267-4150-a520-444f5fc2f27e
md"""
### 2D Example
"""

# ╔═╡ 03d28001-8805-4ba9-8aa3-1f38fa06c219
g2=HelmholtzProject.grid2d_circle_in_square(nref=2)

# ╔═╡ f84999b9-c698-45de-9540-5155a0f1fa4c
gridplot(g2,resolution=(300,300))	

# ╔═╡ 89190b87-80c5-49a1-b341-1c426aabdfdc
function bc2(i,xy) 
    i==2 && return  Dirichlet(),1
    i==4 && return  Dirichlet(),0
    return 0,0
end

# ╔═╡ 0c78c3b2-cbea-4bef-a236-5d94048c963a
function k2(i,xy)
	i==1 && return 1
	i==2 && return 10000
end

# ╔═╡ 7d269042-f893-4c63-be8b-3f241cd8c52e
A2,f2=assemble_reaction_diffusion(g2,b=bc2,k=k2)

# ╔═╡ 6de8aac9-fd79-4bcd-a5a5-0076d6dcb0a9
sol2=A2\f2

# ╔═╡ 805178e7-998f-4929-b455-f81d90e7b760
scalarplot(g2,sol2,resolution=(300,300),levels=7)

# ╔═╡ fa511aa2-52f2-41b1-9ac5-847ce30b37d9
md"""
### 3D Example
"""

# ╔═╡ 77ee9cdb-4519-4750-8657-1bc836f8b81d
X3=collect(0.0:0.05:1)

# ╔═╡ 4ddf06e0-d8b2-4cc9-bfb9-3d465d163995
g3=simplexgrid(X3,X3,X3)

# ╔═╡ fa80f8ff-19b5-41f4-92f0-b465948d65bc
gridplot(g3,resolution=(300,300),zplane=0.0,outlinealpha=0.4)

# ╔═╡ 3f8c9b45-73f3-4887-a578-9d061610226f
function f3(i,xyz)
	return xyz[1]*xyz[2]*xyz[3]
end

# ╔═╡ 4179c8f7-1730-442e-bd56-38153204aaf6
A3,b3=assemble_reaction_diffusion(g3,f=f3)

# ╔═╡ 164eca31-9b98-46ba-8921-ddbbf805f469
sol3=A3\b3

# ╔═╡ 7d06bf0b-62a2-44e2-8ec8-94d4779b6e06
	scalarplot(g3,sol3,resolution=(300,300))

# ╔═╡ Cell order:
# ╠═f9a43662-326a-46b4-a638-39498fb59ce6
# ╟─0f649cee-b670-4c66-9c6f-72e6a93f52a0
# ╟─9eaaa834-0b7b-46a6-a300-a70a5c7af9db
# ╠═f514b55e-92b6-418b-8a29-cf0c425453c8
# ╠═b1e6e0d8-c22d-453d-9deb-cbde27656689
# ╠═1d888171-2b58-4370-8bf2-63905a3d37f9
# ╠═4de4e9c5-ee40-40be-9bed-9f538d45674e
# ╠═9888e062-57f7-434a-9176-0a3cd8230448
# ╠═37f9ca73-d1e6-4afc-b2a3-858a80faa2e5
# ╟─f65baba2-6267-4150-a520-444f5fc2f27e
# ╠═03d28001-8805-4ba9-8aa3-1f38fa06c219
# ╠═f84999b9-c698-45de-9540-5155a0f1fa4c
# ╠═89190b87-80c5-49a1-b341-1c426aabdfdc
# ╠═0c78c3b2-cbea-4bef-a236-5d94048c963a
# ╠═7d269042-f893-4c63-be8b-3f241cd8c52e
# ╠═6de8aac9-fd79-4bcd-a5a5-0076d6dcb0a9
# ╠═805178e7-998f-4929-b455-f81d90e7b760
# ╟─fa511aa2-52f2-41b1-9ac5-847ce30b37d9
# ╠═77ee9cdb-4519-4750-8657-1bc836f8b81d
# ╠═4ddf06e0-d8b2-4cc9-bfb9-3d465d163995
# ╠═fa80f8ff-19b5-41f4-92f0-b465948d65bc
# ╠═3f8c9b45-73f3-4887-a578-9d061610226f
# ╠═4179c8f7-1730-442e-bd56-38153204aaf6
# ╠═164eca31-9b98-46ba-8921-ddbbf805f469
# ╠═7d06bf0b-62a2-44e2-8ec8-94d4779b6e06
