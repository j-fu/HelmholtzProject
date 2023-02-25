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
## 1D case
"""

# ╔═╡ f2d62346-1c3a-41c5-a268-d53874fab699
md"""
``\Omega=(0,1)``, ``-\Delta u - k_0^2 u=\lambda u``, ``u'(0)=u'(1)=0``

``u(x)=cos((k-1)πx)``

``\lambda_k=(k-1)^2\pi^2 - k_0^2``
"""

# ╔═╡ f5a5f39f-4a2b-4891-b893-f4702e826d08
md"""
Finest 1D grid
"""

# ╔═╡ 1a35c2bd-d9b4-4e7b-ab3d-891dcb45267f
maxnref1D = 9

# ╔═╡ 7e4fbf5e-44ed-45bc-86c7-1334d050678f
md"""
Number of eigenvalues to consider
"""

# ╔═╡ 2941e9fd-2c79-427a-b72f-167151569b16
k_0_1d=100+100.0im

# ╔═╡ fdcfe46c-f6e2-4c37-b834-55223f4e664c
nev=10

# ╔═╡ 37bbefdc-c4a9-41cd-bb81-74d8ec3e1b01
eigensolver=jd_eigen

# ╔═╡ 5d93dff6-054b-4a7e-a9fa-6548debd78f5
md"""
Use reference solution from finest grid ?
"""

# ╔═╡ 51c9dfe4-9f07-480f-8447-c8649e86aa56
use_vref_finest=true

# ╔═╡ 608499e5-407b-4924-ba10-97eb2b9bab5b
λref=[(i-1)^2*π^2 - k_0_1d^2 for i=1:nev]

# ╔═╡ e8c56784-bf8a-426b-983e-bd32479bd81a
Xref=0:2.0^-maxnref1D:1

# ╔═╡ b0ee4c1c-5463-4b5d-ad40-3005fce3d1d7
gridref=simplexgrid(Xref)

# ╔═╡ 4af3d30c-6f67-4fcf-b4f6-59d7d1f24f85
Aref,Mref=assemble_helmholtz(gridref,k=k_0_1d,n=(i,x)->1.0,m=mass_matrix_lumped);

# ╔═╡ ae0db901-996e-46e5-b7b9-06cd95a0a2be
vref_finest=[HelmholtzProject.normalize!(map(x->cos((i-1)*π*x),gridref),Mref) for i=1:nev]

# ╔═╡ fefb13a7-a026-48e9-8825-991ebd35d073
# ╠═╡ show_logs = false
begin
	Error_λ = []
	Error_v = []
	h=[]
	for nref=4:maxnref1D
		X=0:2.0^-nref:1
		push!(h,X[2]-X[1])
		grid = simplexgrid(X)
		A,M=assemble_helmholtz(grid,k=k_0_1d,n=(i,x)->1.0,m=mass_matrix_lumped)
		λ,v=eigensolver(A,M,nev=nev,maxiter=10000,tol=1.0e-10); 
		vref=[HelmholtzProject.normalize!(map(x->cos((i-1)*π*x),grid),M) for i=1:nev]
		λ_error=[]
		v_errorl2=[]
		v_errorh1=[]
		for i=1:nev
			push!(λ_error,abs(λ[i] - λref[i]))
			if use_vref_finest
				vinter=HelmholtzProject.interpolate(gridref,v[:,i],grid)
				vl2,vh1 = HelmholtzProject.fenorms(abs.(vinter) - abs.(vref_finest[i]),gridref)
			else
				vl2,vh1 = HelmholtzProject.fenorms(abs.(v[:,i]) - abs.(vref[i]),grid)
			end
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
	error(iλ)=[Error_λ[i][iλ] for i=1:length(Error_λ)]
	for i=2:nev
     	scalarplot!(vis,h,error(i),color=RGB(1-(i-1)/nev,0,(i-1)/nev),clear=false,label="k=$i")
	end
	scalarplot!(vis,h,0.2*h.^2,label="O(h^2)",clear=false,color=:black,linestyle=:dot)
	reveal(vis)
end

# ╔═╡ 744db4dc-2e65-47ce-87e9-b96c6dd0f4aa
let
	vis=GridVisualizer(size=(500,300),xscale=:log,yscale=:log,legend=:lt,title="Eigenvector error")
	error(iλ)=[Error_v[i][iλ] for i=1:length(Error_λ)]
	for i=2:nev
     	scalarplot!(vis,h,error(i),color=RGB(1-(i-1)/nev,0,(i-1)/nev),clear=false,label="k=$i")
	end
	scalarplot!(vis,h,0.2*h.^2,label="O(h^2)",clear=false,color=:black,linestyle=:dot)
	scalarplot!(vis,h,0.05*h,label="O(h)",clear=false,color=:black,linestyle=:dash,linewidth=1)
	reveal(vis)
end

# ╔═╡ 28654c57-e6f6-4fc6-a841-64e9a7c89820
md"""
# 2D Case
"""

# ╔═╡ 4ea0e288-6196-4b39-be00-e0be0ab1c5a3
dirichlet2d=false

# ╔═╡ dc1d2980-0972-4f00-ae97-505ef6681a75
md"""
``\Omega=(0,w_x)\times (0,w_y)``, ``-\Delta u=\lambda u``, ``\partial_n u=0``

``u(x,y)=cos((k-1)πx/w_x)cos((l-1)πy/w_y)``

``\lambda_{kl}=((l-1)^2/w_x^2+(k-1)^2/w_y^2)π^2`` 
"""

# ╔═╡ 2bc9b41b-9cdf-4285-94fc-308ed0542801
wx=1.0; wy=1.25

# ╔═╡ 7ebb814d-e433-492a-b883-255ae414cf37
k_0_2d=10.0+10im

# ╔═╡ 64fd813d-9991-42ba-8f1c-a9f003b663ef
maxnref2d=5

# ╔═╡ b08cd584-2e26-476f-af1c-79f95b2874d6
nev2d=10

# ╔═╡ 1f9e5f6b-6018-4fb3-a697-c913de362288
rect=rect_r

# ╔═╡ 50d7b0e8-af27-416d-9c1b-08f19d0487eb
gridref2d=rect(nref=maxnref2d,w=wx,h=wy)

# ╔═╡ 7f3ae857-be85-44fa-83ee-f600681267f7
gridplot(gridref2d, size=(200,200))

# ╔═╡ 45f52a3d-76ec-4481-bce5-43ac67ab7b85
md"""
We need to sort analytical eigenvalues after their size. 
For this we create a permutation which then is also used to sort the eigenvectors.
"""

# ╔═╡ 8b841921-32bd-4968-9e7a-58058e7f50c7
begin
	off=dirichlet2d ? 1 : 0
 	λref2d=vec([((k+off)^2/wx^2+(l+off)^2/wy^2)*π^2-k_0_2d^2 for k=0:nev2d, l=0:nev2d])
	perm2d=sortperm(abs.(λref2d))
	λref2d=λref2d[perm2d]
end

# ╔═╡ fe090c86-06f2-4c3a-a7c8-eb9b3364d6c2
if dirichlet2d
	b2d(i,x)=Dirichlet(),0
else
	b2d(i,x)=0,0
end

# ╔═╡ e4df6a0f-dd50-43a4-9dcf-36490c43136c
begin
    Aref2d,Mref2d=assemble_helmholtz(gridref2d,k=k_0_2d,
	n=(i,x)->1.0,b=b2d,m=mass_matrix_lumped);
	if dirichlet2d
		Dinv=inv(Diagonal(Aref2d))
		Aref2d=Dinv*Aref2d
		Mref2d=Dinv*Mref2d
	end
end

# ╔═╡ 586d9db6-f145-4bdb-8828-a901d7b34288
eigensolver2=arpack_eigen

# ╔═╡ 996186c6-ac51-454b-85b6-ed5013623fb0
# ╠═╡ show_logs = false
λ2,v2=eigensolver2(Aref2d,Mref2d,nev=nev2d,maxiter=10000,tol=1.0e-14)

# ╔═╡ b2e11d07-e8d6-4e5a-a324-258ba2eae314
iλ2=10

# ╔═╡ a7e3f275-50bf-4847-ade1-2cdea5307d62
fenorms(Aref2d*v2[:,iλ2]-λ2[iλ2]*Mref2d*v2[:,iλ2],gridref2d)

# ╔═╡ 8dc246d7-41a6-493f-9e7d-c2a39826be75
function vref2d(grid,M)
if dirichlet2d
v=vec([HelmholtzProject.normalize!(map((x,y)->sin((k+1)*π*x/wx)*sin((l+1)*π*y/wy),
	grid),M) for k=0:nev2d,l=0:nev2d])
else
v=vec([HelmholtzProject.normalize!(
		map((x,y)->cos(k*π*x/wx)*cos(l*π*y/wy),grid),
		M) for k=0:nev2d,l=0:nev2d])
end
v[perm2d]
end;

# ╔═╡ bebcaf6e-7aed-465d-9587-19828a38c6a4
vref2d_finest=vref2d(gridref2d,Mref2d)

# ╔═╡ f6901680-9c46-4bc4-ba9d-342384f2fed2
EL.grid([
	md"fem λ[$(iλ2)] = $(round(λ2[iλ2],digits=5))"  md"exact λ[$(iλ2)] = $(round(λref2d[iλ2],digits=5))";
	scalarplot(gridref2d,abs.(v2[:,iλ2]),size=(200,200),colormap=:hot)	scalarplot(gridref2d,abs.(vref2d_finest[iλ2]),size=(200,200),colormap=:hot)]
)

# ╔═╡ 4d5addad-d16c-4f80-8db3-62fc0bc874ff
fenorms(Aref2d*vref2d_finest[iλ2]-λref2d[iλ2]*Mref2d*vref2d_finest[iλ2],gridref2d)

# ╔═╡ 97ec3cb0-15fb-4090-8d84-d307362d9648
fenorms(Aref2d*vref2d_finest[iλ2]-λ2[iλ2]*Mref2d*vref2d_finest[iλ2],gridref2d)

# ╔═╡ 1031c32c-4d1a-41fb-a7ae-557532c1391c
use_vref2d_finest=true

# ╔═╡ c87f3096-ff2a-4a80-8686-b22924b7ba90
# ╠═╡ show_logs = false
begin
	Error2d_λ = []
	Error2d_v = []
	h2d=[]
	for nref=2:maxnref2d
		grid = rect(;nref,w=wx,h=wy)
		push!(h2d,hminmax(grid)[1])
		A,M=assemble_helmholtz(grid,k=k_0_2d,n=(i,x)->1.0,m=mass_matrix_lumped)
		λ,v=eigensolver2(A,M,nev=nev2d,maxiter=10000,tol=1.0e-10); 
		vref=vref2d(grid,M)
		λ_error=[]
		v_errorl2=[]
		v_errorh1=[]
		for i=1:nev
			push!(λ_error,abs(λ[i] - λref2d[i]))
			if use_vref2d_finest
				vinter=HelmholtzProject.interpolate(gridref2d,v[:,i],grid,eps=1.0e-10)
				vl2,vh1 = HelmholtzProject.fenorms(abs.(vinter) - abs.(vref2d_finest[i]),gridref2d)
				@info vl2
			else
				vl2,vh1 = HelmholtzProject.fenorms(abs.(v[:,i]) - abs.(vref[i]),grid)
			end
	        push!(v_errorl2,vl2)	
	        push!(v_errorh1,vh1)	
		end
		push!(Error2d_λ, λ_error)
		push!(Error2d_v, v_errorl2)
	end
end

# ╔═╡ 5f99ac9d-37a6-4487-a28d-ece7bc4e9b4d
Error2d_λ

# ╔═╡ 77366ca0-5c30-496d-92a2-1c28445dc07e
let
	vis=GridVisualizer(size=(500,300),xscale=:log,yscale=:log,legend=:lt,title="Eigenvalue error",xlabel="h",ylabel="|λ-λ_ref|")
	error(iλ)=[Error2d_λ[i][iλ] for i=1:length(Error2d_λ)]
	for i=2:nev2d
     	scalarplot!(vis,h2d,error(i),color=RGB(1-(i-1)/nev2d,0,(i-1)/nev2d),clear=false,label="k=$i")
	end
	scalarplot!(vis,h2d,0.2*h2d.^2,label="O(h^2)",clear=false,color=:black,linestyle=:dot)
	reveal(vis)
end

# ╔═╡ 5a8428af-4312-42b2-ad05-040f0e03f544
let
	vis=GridVisualizer(size=(500,300),xscale=:log,yscale=:log,legend=:lt,title="Eigenvector error",xlabel="h",ylabel="|v-v_ref|")
	error(iλ)=[Error2d_v[i][iλ] for i=1:length(Error2d_v)]
	for i=2:nev2d
     	scalarplot!(vis,h2d,error(i),color=RGB(1-(i-1)/nev2d,0,(i-1)/nev2d),clear=false,label="k=$i")
	end
	scalarplot!(vis,h2d,0.8*h2d,label="O(h)",clear=false,color=:black,linestyle=:dash)
scalarplot!(vis,h2d,0.2*h2d.^2,label="O(h^2)",clear=false,color=:black,linestyle=:dot)
		reveal(vis)
end

# ╔═╡ Cell order:
# ╠═8d3aecfa-dcdf-11ec-1933-192af76ee162
# ╟─0f1b4112-3a36-46ee-8bc3-61a1890babdf
# ╟─f2d62346-1c3a-41c5-a268-d53874fab699
# ╟─f5a5f39f-4a2b-4891-b893-f4702e826d08
# ╠═1a35c2bd-d9b4-4e7b-ab3d-891dcb45267f
# ╟─7e4fbf5e-44ed-45bc-86c7-1334d050678f
# ╠═2941e9fd-2c79-427a-b72f-167151569b16
# ╠═fdcfe46c-f6e2-4c37-b834-55223f4e664c
# ╠═37bbefdc-c4a9-41cd-bb81-74d8ec3e1b01
# ╟─5d93dff6-054b-4a7e-a9fa-6548debd78f5
# ╠═51c9dfe4-9f07-480f-8447-c8649e86aa56
# ╠═608499e5-407b-4924-ba10-97eb2b9bab5b
# ╠═e8c56784-bf8a-426b-983e-bd32479bd81a
# ╠═b0ee4c1c-5463-4b5d-ad40-3005fce3d1d7
# ╠═4af3d30c-6f67-4fcf-b4f6-59d7d1f24f85
# ╠═ae0db901-996e-46e5-b7b9-06cd95a0a2be
# ╠═fefb13a7-a026-48e9-8825-991ebd35d073
# ╠═3b55519f-d1ba-4d34-9e1c-3cefb98407e6
# ╠═293748e1-49c7-44e1-99d3-6a0c00531110
# ╠═a7a3911d-f766-41cf-93fb-73871bf344cd
# ╠═744db4dc-2e65-47ce-87e9-b96c6dd0f4aa
# ╟─28654c57-e6f6-4fc6-a841-64e9a7c89820
# ╠═4ea0e288-6196-4b39-be00-e0be0ab1c5a3
# ╟─dc1d2980-0972-4f00-ae97-505ef6681a75
# ╠═2bc9b41b-9cdf-4285-94fc-308ed0542801
# ╠═7ebb814d-e433-492a-b883-255ae414cf37
# ╠═64fd813d-9991-42ba-8f1c-a9f003b663ef
# ╠═b08cd584-2e26-476f-af1c-79f95b2874d6
# ╠═1f9e5f6b-6018-4fb3-a697-c913de362288
# ╠═50d7b0e8-af27-416d-9c1b-08f19d0487eb
# ╠═7f3ae857-be85-44fa-83ee-f600681267f7
# ╟─45f52a3d-76ec-4481-bce5-43ac67ab7b85
# ╠═8b841921-32bd-4968-9e7a-58058e7f50c7
# ╠═fe090c86-06f2-4c3a-a7c8-eb9b3364d6c2
# ╠═e4df6a0f-dd50-43a4-9dcf-36490c43136c
# ╠═586d9db6-f145-4bdb-8828-a901d7b34288
# ╠═996186c6-ac51-454b-85b6-ed5013623fb0
# ╟─f6901680-9c46-4bc4-ba9d-342384f2fed2
# ╠═b2e11d07-e8d6-4e5a-a324-258ba2eae314
# ╠═a7e3f275-50bf-4847-ade1-2cdea5307d62
# ╠═4d5addad-d16c-4f80-8db3-62fc0bc874ff
# ╠═97ec3cb0-15fb-4090-8d84-d307362d9648
# ╠═8dc246d7-41a6-493f-9e7d-c2a39826be75
# ╠═bebcaf6e-7aed-465d-9587-19828a38c6a4
# ╠═1031c32c-4d1a-41fb-a7ae-557532c1391c
# ╠═c87f3096-ff2a-4a80-8686-b22924b7ba90
# ╠═5f99ac9d-37a6-4487-a28d-ece7bc4e9b4d
# ╠═77366ca0-5c30-496d-92a2-1c28445dc07e
# ╠═5a8428af-4312-42b2-ad05-040f0e03f544
