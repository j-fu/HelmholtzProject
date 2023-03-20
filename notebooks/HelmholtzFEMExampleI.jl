### A Pluto.jl notebook ###
# v0.19.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

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
	using LinearAlgebra
	using ExtendableSparse
 	using HelmholtzProject
	default_plotter!(PlutoVista)	
end

# ╔═╡ 3758a738-5fdc-4a42-97ae-095715978b68
md"""
## 1D case
"""

# ╔═╡ f2b05487-5608-4c5f-956a-05ae291f9967
N1=5000

# ╔═╡ 7357d693-5370-4369-bcd2-25e3136f1d45
nev1=20

# ╔═╡ 39e4676c-2da7-4060-bb4f-97212b1e294c
X1=0:1/N1:1

# ╔═╡ d834cc33-31dc-42c3-91f7-1a5349660b98
grid1=simplexgrid(X1)

# ╔═╡ 9716572c-9f11-494d-b34b-4596406213de
n(x)= abs(x-1/2) < 0.1 ? 1 : 0.5 ;

# ╔═╡ b9733cd3-502d-4ff1-912b-a02bfd1fd333
b1(i,coord)= (0, 0) 

# ╔═╡ 1ffd1959-97ee-4c18-8607-93917b28300b
md"""
real(k)= $(@bind rkT Slider(0:0.5:500,show_value=true,default=200))
imag(k)= $(@bind ikT Slider(-50:0.1:50,show_value=true,default=1.0))
"""

# ╔═╡ d1a2d08e-2519-4c2e-b9a8-6a07fd2ad844
kT=rkT+ikT*im;kT^2

# ╔═╡ ddca0556-9b43-4125-989d-e1a0547a3806
md"""
Eigenvalue index $(@bind iT PlutoUI.Slider(1:nev1,show_value=true))
"""

# ╔═╡ 6727d195-5897-40d9-9a66-a6e572736783
viseT1=GridVisualizer(resolution=(600,150),legend=:rt); viseT1

# ╔═╡ 16d201b6-b352-4130-a9f7-4dd702319442
viseT2=GridVisualizer(resolution=(600,150),legend=:rt); viseT2

# ╔═╡ 412b131a-815f-4564-aeb2-61032e0a404a
md"""
## 2D Case
"""

# ╔═╡ 226820e7-d084-4193-88e7-8674c22861e8
N2=101

# ╔═╡ 8d083b2c-4f3e-474a-a154-cc1abd867f6e
X20=collect(range(0,0.4,length=N2÷4))

# ╔═╡ 8c1d6361-b807-4350-a59b-39e9cb88e07f
X21=collect(range(0.4,0.5,length=N2÷4))

# ╔═╡ 69a3041d-09f3-4298-9e35-daf50285ff3f
X2=glue(glue(X20,X21), 1.0.-reverse(glue(X20,X21)))

# ╔═╡ 9a9289f4-6c69-4f6d-95de-9f67a320ef92
circle=true

# ╔═╡ a9253faa-9e95-45ec-bc7c-df750143f737
nref=5

# ╔═╡ 34d212dc-856c-4eb4-a677-ae3858a1ca89
n(x,y)= abs(x-1/2) < 0.1 && abs(y-1/2) < 0.1 ? 1 : 0.5 ;

# ╔═╡ f074a523-68fb-48d0-9623-1066c87508d4
scalarplot(X1,n.(X1),resolution=(600,150))

# ╔═╡ 98aeb541-c547-4235-b373-9941ac6080fc
n1(i,coord)=n(coord[1])

# ╔═╡ 3d8e540d-f5d6-418b-822d-4db372333f84
A1,M1=assemble_helmholtz(grid1,k=kT,n=n1,m=mass_matrix_lumped)

# ╔═╡ 69dc786c-d724-44a2-a803-bd8070b8f243
# ╠═╡ show_logs = false
@time λ1,v1=jd_eigen(A1,M1,nev=nev1,maxiter=10000,tol=1.0e-10);

# ╔═╡ a7f0add2-cfe4-47d5-b3f1-c1563137c455
length(λ1)

# ╔═╡ 52111b28-fc99-413c-a702-88bc252f383e
md"""λ($iT) = $(round(λ1[iT],sigdigits=5))"""

# ╔═╡ 089c2f36-fe58-4450-8aed-cb3f7825b1a2
begin
    @views uet=v1[:,iT]
	l2,h1=HelmholtzProject.fenorms(uet,grid1)
	@info "l2 norm: $(l2), h1 norm: $(h1)"
    scalarplot!(viseT1,
        X1,imag(uet),label="im",color=:blue)
    scalarplot!(viseT1,
        X1,real(uet),label="re",clear=false,color=:red,show=true)
end

# ╔═╡ e3edbbd9-19e4-41ee-8cba-05575b05c1b1
scalarplot!(viseT2,
        X1,abs.(uet),label="abs",color=:green,show=true)

# ╔═╡ 1f02a3da-70c0-487b-b999-7de1a4c88375
if circle
  grid2=HelmholtzProject.grid2d_circle_in_square(nref=nref)
	n2(i,coord)= i==1 ? 0.5 : 1 
else
	grid2=simplexgrid(X2,X2)
	n2(i,coord)=n(coord[1],coord[2])
end

# ╔═╡ caedced1-cf90-475e-88ba-a5880762fda7
gridplot(grid2,size=(300,300))

# ╔═╡ 3923b8c9-65ea-4427-bd36-ed61d380299e
nev=20

# ╔═╡ 2ab796c4-ad9c-469f-a358-6aa48fbcfba0
function plotev(λ0)
	vis=GridVisualizer(size=(500,200),legend=:rt,xlabel="eigenvalue index")
	λ=λ0[1:min(nev,length(λ0))]
	r=1:length(λ)
	scalarplot!(vis,r,abs.(λ),label="|λ|",color=:green)
	scalarplot!(vis,r,real.(λ),label="Re(λ)",color=:red,clear=false)
	scalarplot!(vis,r,imag.(λ),label="Im(λ)",color=:blue,clear=false)
	reveal(vis)
end

# ╔═╡ 5d18de00-a949-49a7-9a9e-650da7c2bdbb
plotev(λ1)

# ╔═╡ b18bc1b1-6e74-4d64-865d-966384163ecf
md"""
real(k)= $(@bind rkT2 confirm(Slider(0:0.5:500,show_value=true,default=200)))
imag(k)= $(@bind ikT2 confirm(Slider(-50:0.1:50,show_value=true,default=1.0)))
"""

# ╔═╡ 0737c4b4-f6fc-4578-86be-0dcfffc37d16
kT2=rkT2+ikT2*im;kT2^2

# ╔═╡ 9c3359c6-afa4-4665-9933-98b8ebcf0bc9
A2,M2=assemble_helmholtz(grid2,k=kT2,n=n2)

# ╔═╡ 5600a195-db61-4bf7-815d-1392dcea771f
# ╠═╡ show_logs = false
@time λ,v=jd_eigen(A2,M2;nev=nev,maxiter=500);

# ╔═╡ ba542ed8-5408-4b47-9877-0d6651aa472c
md"""
Eigenvalue index $(@bind iT2 PlutoUI.Slider(1:length(λ),show_value=true))
"""

# ╔═╡ 5aec611b-94a9-420b-a71c-64b09ee7e373
λ[iT2]

# ╔═╡ 9dee97de-c7c8-4209-a2bd-c22a3e394b67
vis=GridVisualizer(dim=2,size=(500,300))

# ╔═╡ db8521f3-9ceb-4f63-8032-9da454d54577
let
	@info iT2
	scalarplot!(vis,grid2,abs.(v[:,iT2]),label="abs",colormap=:hot,show=true,levels=0)
end

# ╔═╡ 2e873ab7-08e1-4d8b-acbe-dcaec4218ba9
plotev(λ)

# ╔═╡ 2a8bf655-66da-4567-9800-ad8f6826d970
function bary!(λ,invA,L2G,x)
    mapderiv!(invA,L2G,nothing)
    fill!(λ ,0)
    for j = 1 : length(x)
        dj=x[j] - L2G.b[j]
        for  k = 1 : length(x)
            λ[k] += invA[j,k] * dj
        end
    end
#    ExtendableGrids.postprocess_xreftest!(λ,Triangle2D)
end


# ╔═╡ a082642c-4809-42fe-a912-3fb60646ff8d
g0=simplexgrid(0:0.1:1,0:0.1:1)

# ╔═╡ dbf1cdb7-a35a-463a-8c7d-77986a986a5f
u0=map((x,y)->x+2y,g0)

# ╔═╡ 5eebbd19-913f-48cf-a644-92bd3e88ac1d
gx=simplexgrid(0:0.05:1,0:0.05:1)

# ╔═╡ af0d74c4-06f8-4867-ac58-c799c9a0fca8
function grid2d(;nref=0)
b=SimplexGridBuilder(Generator=Triangulate)
	p1=point!(b,0,0)
	p2=point!(b,1,0)
	p3=point!(b,1,1)
	p4=point!(b,0,1)

	# Specify outer boundary
	facetregion!(b,1)
	facet!(b,p1,p2)
	facetregion!(b,2)
	facet!(b,p2,p3)
	facetregion!(b,3)
	facet!(b,p3,p4)	
	facetregion!(b,4)
	facet!(b,p4,p1)	
	
	simplexgrid(b,maxvolume=4.0^(-nref))
end


# ╔═╡ 2ad3ba71-be38-41de-b778-0761b56fb75a
g1=grid2d(;nref=4)

# ╔═╡ 46f8d6f0-3fb5-40bf-aacf-498e17e3f6fa
u1=HelmholtzProject.interpolate(g1,u0,g0)

# ╔═╡ db1151e3-c65b-4442-8f7c-1959278d805d
scalarplot(g1,u1,size=(500,300))

# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─3758a738-5fdc-4a42-97ae-095715978b68
# ╠═f2b05487-5608-4c5f-956a-05ae291f9967
# ╠═7357d693-5370-4369-bcd2-25e3136f1d45
# ╠═39e4676c-2da7-4060-bb4f-97212b1e294c
# ╠═d834cc33-31dc-42c3-91f7-1a5349660b98
# ╠═9716572c-9f11-494d-b34b-4596406213de
# ╠═f074a523-68fb-48d0-9623-1066c87508d4
# ╠═98aeb541-c547-4235-b373-9941ac6080fc
# ╠═b9733cd3-502d-4ff1-912b-a02bfd1fd333
# ╠═3d8e540d-f5d6-418b-822d-4db372333f84
# ╠═69dc786c-d724-44a2-a803-bd8070b8f243
# ╠═a7f0add2-cfe4-47d5-b3f1-c1563137c455
# ╟─1ffd1959-97ee-4c18-8607-93917b28300b
# ╠═d1a2d08e-2519-4c2e-b9a8-6a07fd2ad844
# ╟─ddca0556-9b43-4125-989d-e1a0547a3806
# ╟─52111b28-fc99-413c-a702-88bc252f383e
# ╠═6727d195-5897-40d9-9a66-a6e572736783
# ╟─16d201b6-b352-4130-a9f7-4dd702319442
# ╠═089c2f36-fe58-4450-8aed-cb3f7825b1a2
# ╟─e3edbbd9-19e4-41ee-8cba-05575b05c1b1
# ╟─5d18de00-a949-49a7-9a9e-650da7c2bdbb
# ╟─2ab796c4-ad9c-469f-a358-6aa48fbcfba0
# ╟─412b131a-815f-4564-aeb2-61032e0a404a
# ╠═226820e7-d084-4193-88e7-8674c22861e8
# ╠═8d083b2c-4f3e-474a-a154-cc1abd867f6e
# ╠═8c1d6361-b807-4350-a59b-39e9cb88e07f
# ╠═69a3041d-09f3-4298-9e35-daf50285ff3f
# ╠═9a9289f4-6c69-4f6d-95de-9f67a320ef92
# ╠═a9253faa-9e95-45ec-bc7c-df750143f737
# ╠═1f02a3da-70c0-487b-b999-7de1a4c88375
# ╠═caedced1-cf90-475e-88ba-a5880762fda7
# ╠═34d212dc-856c-4eb4-a677-ae3858a1ca89
# ╠═9c3359c6-afa4-4665-9933-98b8ebcf0bc9
# ╠═0737c4b4-f6fc-4578-86be-0dcfffc37d16
# ╠═3923b8c9-65ea-4427-bd36-ed61d380299e
# ╟─b18bc1b1-6e74-4d64-865d-966384163ecf
# ╠═5600a195-db61-4bf7-815d-1392dcea771f
# ╟─ba542ed8-5408-4b47-9877-0d6651aa472c
# ╠═5aec611b-94a9-420b-a71c-64b09ee7e373
# ╠═9dee97de-c7c8-4209-a2bd-c22a3e394b67
# ╠═db8521f3-9ceb-4f63-8032-9da454d54577
# ╟─2e873ab7-08e1-4d8b-acbe-dcaec4218ba9
# ╠═2a8bf655-66da-4567-9800-ad8f6826d970
# ╠═a082642c-4809-42fe-a912-3fb60646ff8d
# ╠═dbf1cdb7-a35a-463a-8c7d-77986a986a5f
# ╠═5eebbd19-913f-48cf-a644-92bd3e88ac1d
# ╠═af0d74c4-06f8-4867-ac58-c799c9a0fca8
# ╠═2ad3ba71-be38-41de-b778-0761b56fb75a
# ╠═46f8d6f0-3fb5-40bf-aacf-498e17e3f6fa
# ╠═db1151e3-c65b-4442-8f7c-1959278d805d
