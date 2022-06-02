### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
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
	default_plotter!(PlutoVista)	

end

# ╔═╡ 39e4676c-2da7-4060-bb4f-97212b1e294c
X1=0:0.01:1

# ╔═╡ d834cc33-31dc-42c3-91f7-1a5349660b98
grid1=simplexgrid(X1)

# ╔═╡ 98aeb541-c547-4235-b373-9941ac6080fc
r1(i,coord)=100.0

# ╔═╡ b9733cd3-502d-4ff1-912b-a02bfd1fd333
b1(i,coord)= i==1 ? (Dirichlet(), 1) : (Dirichlet(), 0) 

# ╔═╡ 78a226ef-81e7-4ebc-af15-31099ede6059
A1,f1=assemble_reaction_diffusion(grid1,r=r1,b=b1)

# ╔═╡ 24076386-b486-4aab-990a-d5f4d57c088d
sol1=A1\f1

# ╔═╡ da7f2d46-e59e-4a9b-a7c5-6719cee4d657
scalarplot(grid1,sol1,resolution=(600,300),color=:red)

# ╔═╡ 6430f557-62b4-4ae3-a7c6-4b8eb7fc73d3
grid2=HelmholtzProject.grid2d_circle_in_square(nref=2)

# ╔═╡ 51fca0bd-97aa-4767-863e-1f287533e3cc
function b2(i,coord)
	i==1 && return Dirichlet(),0
	i==3 && return Dirichlet(),1
	0,0
end

# ╔═╡ bed090c5-e5ec-45db-9855-f314a018b628
function r2(i,coord)
	i==2 && return 1000
	return 0
end

# ╔═╡ 15871888-da43-453a-9dbe-748fb97eaded
A2,f2=assemble_reaction_diffusion(grid2,b=b2,r=r2)

# ╔═╡ 2ae80cec-16cd-4312-97d4-9f37f05d342c
sol2=A2\f2

# ╔═╡ 41ae604a-f138-4e39-840b-39d74b61278e
scalarplot(grid2,sol2,resolution=(300,300))

# ╔═╡ 9529497b-f4de-4942-9471-86ffde4d91ba
md"""
## Convergence test
"""

# ╔═╡ b69f9c93-fe7c-4079-8c59-8707bbeda3ad
k=0;l=1;r=0.1

# ╔═╡ 8c6d4732-db37-41fe-a1d7-899c05a35ec4
rx(i,coord)=r

# ╔═╡ 7dc3857f-a303-4ab1-9026-e02ca3bda059
bx(i,coord)=0.0,0.0

# ╔═╡ 32c58e8d-5da5-41d7-8e87-c51b71c49c33
u0(x,y)=cos(k*π*x)*cos(l*π*y)

# ╔═╡ 4f6dad33-5bcd-4616-a4c8-a25801ceaf55
fx(i,coord)=((k^2+l^2)*π^2+r)*u0(coord...)

# ╔═╡ cee15f5a-44af-45ed-93fc-2eb77189dc0e
maxref=7

# ╔═╡ 12cf1836-1e68-47c0-94bc-f086a8158e67
begin
	l2=[]
	h1=[]
	h=[]
	glast=nothing
	ulast=nothing
	for nref=1:maxref
	    gridt=unitsquare_r(;nref)
  	 	At,ft=assemble_reaction_diffusion(gridt,r=rx,b=bx,f=fx,m=mass_matrix_lumped)
		ut=At\ft
		el2,eh1=fenorms(ut-map(u0,gridt),gridt)
		push!(l2,el2)
		push!(h1,eh1)
		push!(h,hminmax(gridt)[1])
		global glast=gridt
		global ulast=ut
	end
end

# ╔═╡ a5710247-6a10-4cdc-a928-bc73885169e5
let
	vis=GridVisualizer(size=(500,300),xscale=:log,yscale=:log,legend=:lt,xlabel="h",ylabel="error")
	scalarplot!(vis,h,l2,label="l2",color=:green)
	scalarplot!(vis,h,h.^2,label="O(h^2)",color=:green,clear=false,linestyle=:dot)
	scalarplot!(vis,h,h1,label="h1",color=:red,clear=false,linestyle=:none)
	scalarplot!(vis,h,h,label="O(h)",color=:red,clear=false,linestyle=:dot)
	reveal(vis)
end

# ╔═╡ 2b0b57ff-50db-4619-9579-df34e4f0adb6
scalarplot(glast,ulast, colormap=:bwr,size=(300,300))

# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╠═39e4676c-2da7-4060-bb4f-97212b1e294c
# ╠═d834cc33-31dc-42c3-91f7-1a5349660b98
# ╠═98aeb541-c547-4235-b373-9941ac6080fc
# ╠═b9733cd3-502d-4ff1-912b-a02bfd1fd333
# ╠═78a226ef-81e7-4ebc-af15-31099ede6059
# ╠═24076386-b486-4aab-990a-d5f4d57c088d
# ╠═da7f2d46-e59e-4a9b-a7c5-6719cee4d657
# ╠═6430f557-62b4-4ae3-a7c6-4b8eb7fc73d3
# ╠═51fca0bd-97aa-4767-863e-1f287533e3cc
# ╠═bed090c5-e5ec-45db-9855-f314a018b628
# ╠═15871888-da43-453a-9dbe-748fb97eaded
# ╠═2ae80cec-16cd-4312-97d4-9f37f05d342c
# ╠═41ae604a-f138-4e39-840b-39d74b61278e
# ╠═9529497b-f4de-4942-9471-86ffde4d91ba
# ╠═b69f9c93-fe7c-4079-8c59-8707bbeda3ad
# ╠═8c6d4732-db37-41fe-a1d7-899c05a35ec4
# ╠═7dc3857f-a303-4ab1-9026-e02ca3bda059
# ╠═32c58e8d-5da5-41d7-8e87-c51b71c49c33
# ╠═4f6dad33-5bcd-4616-a4c8-a25801ceaf55
# ╠═cee15f5a-44af-45ed-93fc-2eb77189dc0e
# ╠═12cf1836-1e68-47c0-94bc-f086a8158e67
# ╠═a5710247-6a10-4cdc-a928-bc73885169e5
# ╠═2b0b57ff-50db-4619-9579-df34e4f0adb6
