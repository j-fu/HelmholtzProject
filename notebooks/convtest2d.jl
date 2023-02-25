### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 8d3aecfa-dcdf-11ec-1933-192af76ee162
begin 
    using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise
	using ExtendableGrids
	using PlutoVista
	using PlutoUI
	using GridVisualize
	using SimplexGridFactory
	using Triangulate
	using ExtendableSparse
 	using HelmholtzProject
	using LinearAlgebra
	using ForwardDiff
	using Colors
	
	EL=PlutoUI.ExperimentalLayout
	default_plotter!(PlutoVista)	
end

# ╔═╡ 0f1b4112-3a36-46ee-8bc3-61a1890babdf
md"""
## 2D case
"""

# ╔═╡ ab26b33c-6507-4c1c-b1be-0d3d01d795b1
xscale=yscale=1.0

# ╔═╡ 9cb51d12-5117-426a-947d-0db0a11c7f77
gridplot(HelmholtzProject.grid2d_preprint_triangulate(;nref=1))

# ╔═╡ a4b01d8c-bae3-4798-b192-ecb53fb473da
function n(ireg,coord)
	if ireg==1
		return 2.0
	elseif ireg==2
		return 3.0
	else
		return 1.0
	end
end

# ╔═╡ f5a5f39f-4a2b-4891-b893-f4702e826d08
md"""
Finest 1D grid
"""

# ╔═╡ 1a35c2bd-d9b4-4e7b-ab3d-891dcb45267f
maxnref=3

# ╔═╡ b9a997db-dba6-403b-b6be-c5f0ae593270
md"""
Eigenvalue shift
"""

# ╔═╡ 76d29ce0-0bb1-415a-94a3-556fb9af0054
k_0=2*π/(1.55)

# ╔═╡ 7e4fbf5e-44ed-45bc-86c7-1334d050678f
md"""
Number of eigenvalues to consider
"""

# ╔═╡ fdcfe46c-f6e2-4c37-b834-55223f4e664c
nev=4

# ╔═╡ 37bbefdc-c4a9-41cd-bb81-74d8ec3e1b01
eigensolver=arpack_eigen

# ╔═╡ c874521d-ea18-4dd5-91d9-39e90f9383b2
md"""
Calculate reference solution with arpack, as we figured before that  this converges better
"""

# ╔═╡ b0ee4c1c-5463-4b5d-ad40-3005fce3d1d7
gridref=HelmholtzProject.grid2d_preprint_triangulate(;nref=maxnref)

# ╔═╡ 4af3d30c-6f67-4fcf-b4f6-59d7d1f24f85
Aref,Mref=assemble_helmholtz(gridref;k=k_0,n,m=mass_matrix_lumped);

# ╔═╡ da3d5673-585b-4362-b34f-c85a99957d1c
λref,vref=arpack_eigen(Aref,Mref,nev=nev,maxiter=10000,tol=1.0e-7); 

# ╔═╡ 608499e5-407b-4924-ba10-97eb2b9bab5b
λref

# ╔═╡ fefb13a7-a026-48e9-8825-991ebd35d073
# ╠═╡ show_logs = false
begin
	Error_λ = []
	Error_v = []
	h=[]
	for nref=0:maxnref-2
		@info "nref=$(nref)"
		grid = HelmholtzProject.grid2d_preprint_triangulate(;nref=nref,xscale,yscale)
		push!(h,hminmax(grid)[1])
		A,M=assemble_helmholtz(grid;k=k_0,n,m=mass_matrix_lumped)
		λ,v=eigensolver(A,M,nev=nev,maxiter=10000,tol=1.0e-10); 
		λ_error=[]
		v_errorl2=[]
		v_errorh1=[]
		for i=1:nev
			push!(λ_error,abs(λ[i] - λref[i]))
			vinter=HelmholtzProject.interpolate(gridref,v[:,i],grid)
			vl2,vh1 = HelmholtzProject.fenorms(abs.(vinter) - abs.(vref[:,i]),gridref)
	        push!(v_errorl2,vl2)	
	        push!(v_errorh1,vh1)	
		end
		push!(Error_λ, λ_error)
		push!(Error_v, v_errorl2)
	end
end

# ╔═╡ 3b55519f-d1ba-4d34-9e1c-3cefb98407e6
Error_λ

# ╔═╡ 293748e1-49c7-44e1-99d3-6a0c00531110
Error_v

# ╔═╡ a7a3911d-f766-41cf-93fb-73871bf344cd
let
	vis=GridVisualizer(size=(500,300),xscale=:log,yscale=:log,legend=:lt,title="Eigenvalue error",xlabel="h",ylabel="|λ-λ_ref|")
	error(iλ)=[Error_λ[i][iλ]/abs(λref[iλ]) for i=1:length(Error_λ)]
	for i=1:nev
     	scalarplot!(vis,h,error(i),color=RGB(1-(i-1)/nev,0,(i-1)/nev),clear=false,label="k=$i")
	end
	scalarplot!(vis,h,0.4*h.^1.5,label="O(h^(3/2))",clear=false,color=:black,linestyle=:dot)
	reveal(vis)
end

# ╔═╡ 744db4dc-2e65-47ce-87e9-b96c6dd0f4aa
let
	vis=GridVisualizer(size=(500,300),xscale=:log,yscale=:log,legend=:lt,title="Eigenvector error")
	error(iλ)=[Error_v[i][iλ] for i=1:length(Error_λ)]
	for i=1:nev
     	scalarplot!(vis,h,error(i),color=RGB(1-(i-1)/nev,0,(i-1)/nev),clear=false,label="k=$i")
	end
	scalarplot!(vis,h,0.2*h.^2,label="O(h^2)",clear=false,color=:black,linestyle=:dot)
	scalarplot!(vis,h,0.05*h,label="O(h)",clear=false,color=:black,linestyle=:dash,linewidth=1)
	reveal(vis)
end

# ╔═╡ 8645f6a9-6b33-4376-b98f-94066b8282d5
EL.grid([
	scalarplot(gridref,abs.(vref[:,1]),colormap=:hot,size=(300,300)) scalarplot(gridref,abs.(vref[:,2]),colormap=:hot,size=(300,300))
		scalarplot(gridref,abs.(vref[:,3]),colormap=:hot,size=(300,300)) scalarplot(gridref,abs.(vref[:,4]),colormap=:hot,size=(300,300))
])

# ╔═╡ Cell order:
# ╠═8d3aecfa-dcdf-11ec-1933-192af76ee162
# ╠═0f1b4112-3a36-46ee-8bc3-61a1890babdf
# ╠═ab26b33c-6507-4c1c-b1be-0d3d01d795b1
# ╠═9cb51d12-5117-426a-947d-0db0a11c7f77
# ╠═a4b01d8c-bae3-4798-b192-ecb53fb473da
# ╟─f5a5f39f-4a2b-4891-b893-f4702e826d08
# ╠═1a35c2bd-d9b4-4e7b-ab3d-891dcb45267f
# ╟─b9a997db-dba6-403b-b6be-c5f0ae593270
# ╠═76d29ce0-0bb1-415a-94a3-556fb9af0054
# ╟─7e4fbf5e-44ed-45bc-86c7-1334d050678f
# ╠═fdcfe46c-f6e2-4c37-b834-55223f4e664c
# ╠═37bbefdc-c4a9-41cd-bb81-74d8ec3e1b01
# ╠═608499e5-407b-4924-ba10-97eb2b9bab5b
# ╟─c874521d-ea18-4dd5-91d9-39e90f9383b2
# ╠═b0ee4c1c-5463-4b5d-ad40-3005fce3d1d7
# ╠═4af3d30c-6f67-4fcf-b4f6-59d7d1f24f85
# ╠═da3d5673-585b-4362-b34f-c85a99957d1c
# ╠═fefb13a7-a026-48e9-8825-991ebd35d073
# ╠═3b55519f-d1ba-4d34-9e1c-3cefb98407e6
# ╠═293748e1-49c7-44e1-99d3-6a0c00531110
# ╠═a7a3911d-f766-41cf-93fb-73871bf344cd
# ╠═744db4dc-2e65-47ce-87e9-b96c6dd0f4aa
# ╠═8645f6a9-6b33-4376-b98f-94066b8282d5
