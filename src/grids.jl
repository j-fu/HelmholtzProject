### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ ba562efe-1157-11ec-3ed6-97e95cd04353
if isdefined(Main,:PlutoRunner)
    # Ensure not using `,` as floating point decimal delimiter
    # in certain language enviromnents.
    ENV["LC_NUMERIC"]="C"
	ENV["MPLBACKEND"]="agg"
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
    using Revise
	using PlutoUI
	using SimplexGridFactory,GridVisualize,ExtendableGrids,Triangulate,LinearAlgebra
	import PyPlot
	PyPlot.svg(true)
	GridVisualize.default_plotter!(PyPlot)
	inpluto=true
else
	inpluto=false
end;

# ╔═╡ 582b68c7-bf18-47eb-a083-f58d2bdd60ab
md"""
The contents of this pluto notebook is made available as part of the HelmholtzProject project-package. In that case all pluto specific actions and environment activations are disabled. 
"""

# ╔═╡ 9ca39b3c-fca9-4a53-adea-61caa83c05fd
md""" By default, this notebook uses PyPlot as plotting backend for GridVisualize."""

# ╔═╡ 8e564f00-735d-4fca-821d-50fce8e84d79
if inpluto TableOfContents()  end

# ╔═╡ e2ca7b02-4233-4e72-ab26-0f0dcb7bd601
md"""
## Circle in square
"""

# ╔═╡ 6c0000d3-1ef8-43b4-86d6-e65a63577ae3
function grid2d_circle_in_square(;nref=0)
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


	cellregion!(b,2)
        maxvolume!(b,0.05*4.0^(-nref))
	regionpoint!(b,0.5,0.5)
    
	cellregion!(b,1)
        maxvolume!(b,4.0^(-nref))
	regionpoint!(b,0.9,0.9)
        facetregion!(b,5)
	circle!(b,(0.5,0.5),0.25,n=10*2^nref)

	simplexgrid(b)
end


# ╔═╡ 408a3da2-0751-4e78-a8a3-d235f76ea13f
if inpluto  gridplot(grid2d_circle_in_square(;nref=2),resolution=(300,300)) end

# ╔═╡ 87fa7938-4ee7-4749-9538-1116387b808f
md"""##  Grid (preprint)"""

# ╔═╡ 3b62a531-b4a2-4460-926d-03d869649139
md"""
### Default Scaling parameters
"""

# ╔═╡ 73a8bedc-09ee-4d81-adb4-6c02ca1ad86b
begin
	const xscale = 1
	const yscale = 1
	const maxvol = 0.3
	const minvol = 0.01
	const hmin=0.1 * xscale
	const hmax=0.3
end

# ╔═╡ 32099f08-5181-4f60-b9cc-88fed3201a33
md"""
### Unstructured grid using Triangulate
"""

# ╔═╡ 8487f01e-6114-40c3-9d8a-92b3676d4c4c
function builder2d_preprint_triangulate(;xscale=xscale,
		yscale=yscale,maxvol=maxvol,minvol=minvol,nref=0,hmin=hmin,hmax=hmax)
	b=SimplexGridBuilder(Generator=Triangulate;tol=1.0e-10)

	#  Specify points
	p1=point!(b,xscale*0.0,yscale*0.0)
	p2=point!(b,xscale*3.2,yscale*0.0)
	p3=point!(b,xscale*3.2,yscale*3.2)
	p4=point!(b,xscale*0.0,yscale*3.2)

	# Specify outer boundary
	facetregion!(b,1)
	facet!(b,p1,p2)
	facetregion!(b,2)
	facet!(b,p2,p3)
	facetregion!(b,3)
	facet!(b,p3,p4)
	facetregion!(b,4)
	facet!(b,p1,p4)

	# Activate unsuitable callback
	#options!(b,unsuitable=unsuitable)

	# Specify interior boundary
	#GREEN
	facetregion!(b,5)
	facet!(b,point!(b,xscale*0.0,yscale*1.6),point!(b,xscale*0.8,yscale*1.6))
	facet!(b,point!(b,xscale*0.8,yscale*1.6),point!(b,xscale*0.8,yscale*2.4))
	facet!(b,point!(b,xscale*0.8,yscale*2.4),point!(b,xscale*1.5,yscale*2.4))
	facet!(b,point!(b,xscale*1.5,yscale*2.4),point!(b,xscale*1.5,yscale*1.6))
	facet!(b,point!(b,xscale*1.5,yscale*1.6),point!(b,xscale*1.7,yscale*1.6))
	facet!(b,point!(b,xscale*1.7,yscale*1.6),point!(b,xscale*1.7,yscale*2.4))
	facet!(b,point!(b,xscale*1.7,yscale*2.4),point!(b,xscale*2.4,yscale*2.4))
	facet!(b,point!(b,xscale*2.4,yscale*2.4),point!(b,xscale*2.4,yscale*1.6))
	facet!(b,point!(b,xscale*2.4,yscale*1.6),point!(b,xscale*3.2,yscale*1.6))
	#artificial interior boundaries for rectangular shape
	facet!(b,point!(b,xscale*0.8,yscale*1.6),point!(b,xscale*1.5,yscale*1.6))
	facet!(b,point!(b,xscale*1.7,yscale*1.6),point!(b,xscale*2.4,yscale*1.6))

	#BLUE
	facetregion!(b,6)
	facet!(b,point!(b,xscale*0.0,yscale*1.2),point!(b,xscale*3.2,yscale*1.2))
	facet!(b,point!(b,xscale*0.0,yscale*1.4),point!(b,xscale*3.2,yscale*1.4))

	# specify regions
	#substrate
	cellregion!(b,1)
	maxvolume!(b,0.1*4.0^(-nref)*xscale^2)
	regionpoint!(b,xscale*1.0,yscale*2.0)
	regionpoint!(b,xscale*1.0,yscale*1.0)
	regionpoint!(b,xscale*2.0,yscale*2.0)
	regionpoint!(b,xscale*2.0,yscale*1.5)

	# cladding waveguide
	cellregion!(b,2)
	maxvolume!(b,minvol*4.0^(-nref)*xscale^2)
	regionpoint!(b,xscale*2.5,yscale*1.3)

	#air
	cellregion!(b,3)
	maxvolume!(b,maxvol*4.0^(-nref)*xscale^2)
	regionpoint!(b,xscale*2.0,yscale*3.0)



	# Coarse elements in upper left region #1
	#cellregion!(b,1)
	#maxvolume!(b,0.1)
	#regionpoint!(b,0.1,0.5)

	# Fine elements in lower right region #2
	#cellregion!(b,2)
	#maxvolume!(b,0.01)
	#regionpoint!(b,0.9,0.5)

	b
end;

# ╔═╡ 25c049e7-c9c4-4326-a501-34a16d48df7c
grid2d_preprint_triangulate(;kwargs...)=simplexgrid(builder2d_preprint_triangulate(;kwargs...));

# ╔═╡ 626bbd8b-4818-4ad7-9924-0f274ba55561
if inpluto  builder4=builder2d_preprint_triangulate() end;

# ╔═╡ 0deab116-197a-4cfc-a94a-529adfa4ce4f
if inpluto  builderplot(builder4,Plotter=PyPlot,resolution=(750,700)) end

# ╔═╡ 3375c6fe-0ffc-4a80-9f73-245d3f8a4abd
md"""
Unstructured grid using SimplexGrid.
"""

# ╔═╡ 92b4efbe-0bac-4a47-89e5-e482b3ea96be
if inpluto  gridplot(simplexgrid(builder4), resolution=(400,300),linewidth=0.5) end

# ╔═╡ 2a6e529c-cf62-482b-bc99-f0d705cf0ebd
md"""
### Structured grid using SimplexGrid
"""

# ╔═╡ 3db6ebb1-dc8c-4538-bf23-2c7d0b06ad09
function grid2d_preprint_xy(;xscale=xscale,
		yscale=yscale,maxvol=maxvol,hmin=hmin,hmax=hmax)
	l = Int(ceil(1/hmin))

    XA = range(xscale*0.0,xscale*0.8;length=l)
	XB = range(xscale*0.8,xscale*1.5;length=l)
	XC = range(xscale*1.5,xscale*1.7;length=Int(ceil(l/3)))
	XD = range(xscale*1.7,xscale*2.4;length=l)
	XE = range(xscale*2.4,xscale*3.2;length=l)
	
    X = glue(XA,glue(XB,glue(XC,glue(XD,XE))))
	
    YA = geomspace(yscale*0.0,yscale*1.2,hmax,hmin)
	YB = range(yscale*1.2,yscale*1.4;length=l)
	YC = geomspace(yscale*1.4,yscale*1.6,hmin,hmin)
	YD = geomspace(yscale*1.6,yscale*2.4,hmin,hmin)
	YE = geomspace(yscale*2.4,yscale*3.2,hmin,hmax)

	Y = glue(YA,glue(YB,glue(YC,glue(YD,YE))))
	gridnew = simplexgrid(X,Y)

	#unterer Teil
	cellmask!(gridnew, [xscale*0.0,yscale*0.0], [xscale*3.2,YA[end]], 2)

	#wave Teil
	cellmask!(gridnew, [xscale*0.0,YA[end]], [xscale*3.2,YB[end]], 3)

	# Streifen oberhalb wave
	cellmask!(gridnew, [xscale*0.0,YB[end]], [xscale*3.2,YC[end]], 2)
	cellmask!(gridnew, [XA[end],YC[end]],[XB[end],YD[end]], 2)
	cellmask!(gridnew, [XC[end],YC[end]],[XD[end],YD[end]], 2)
    gridnew
end

# ╔═╡ 2878ae3f-8e4b-4f76-a194-d01b8993f6a8
if inpluto
    gridplot(grid2d_preprint_xy(),resolution=(600,400),	linewidth=0.5,legend=:lt)
end

# ╔═╡ b9d2b069-fbd5-4d0a-a4b7-c87cdff5a03d
md"""
Z same achsis as Y, but only using ranges (because still not sure how to handle hmin, hmax)
"""

# ╔═╡ b1c33203-6a8d-4748-87ba-9dfc237f90fc
function grid2d_preprint_xz(;xscale=xscale,
		yscale=yscale,maxvol=maxvol,hmin=hmin,hmax=hmax)

	l = Int(ceil(1/hmin))

    XA = range(xscale*0.0,xscale*0.8;length=l)
	XB = range(xscale*0.8,xscale*1.5;length=l)
	XC = range(xscale*1.5,xscale*1.7;length=Int(ceil(l/3)))
	XD = range(xscale*1.7,xscale*2.4;length=l)
	XE = range(xscale*2.4,xscale*3.2;length=l)

	X = glue(XA,glue(XB,glue(XC,glue(XD,XE))))
	
	ZA = range(yscale*0.0,yscale*1.2;length=l)
	ZB = range(yscale*1.2,yscale*1.4;length=l)
	ZC = range(yscale*1.4,yscale*1.6;length=Int(floor(l/2)))
	ZD = range(yscale*1.6,yscale*2.4,length=l)
	ZE = range(yscale*2.4,yscale*3.2,length=l)
    
	Z = glue(ZA,glue(ZB,glue(ZC,glue(ZD,ZE))))

	gridnewZ = simplexgrid(X,Z)
	
	#unterer Teil
	cellmask!(gridnewZ, [xscale*0.0,yscale*0.0], [xscale*3.2,ZA[end]], 2)

	#wave Teil
	cellmask!(gridnewZ, [xscale*0.0,ZA[end]], [xscale*3.2,ZB[end]], 3)

	# Streifen oberhalb wave
	cellmask!(gridnewZ, [xscale*0.0,ZB[end]], [xscale*3.2,ZC[end]], 2)
	cellmask!(gridnewZ, [XA[end],ZC[end]],[XB[end],ZD[end]], 2)
	cellmask!(gridnewZ, [XC[end],ZC[end]],[XD[end],ZD[end]], 2)

end

# ╔═╡ f58d0c1d-6009-4853-b288-5f1de2d90f40
if inpluto
    gridplot(grid2d_preprint_xz(),resolution=(600,400),
	     linewidth=0.5,legend=:lt)
end

# ╔═╡ 5bd9792b-04eb-4116-a85d-d6201acba562
function rect_t(;nref=0,w=1.0,h=1.0)
b=SimplexGridBuilder(Generator=Triangulate)
	p1=point!(b,0,0)
	p2=point!(b,w,0)
	p3=point!(b,w,h)
	p4=point!(b,0,h)

	# Specify outer boundary
	facetregion!(b,1)
	facet!(b,p1,p2)
	facetregion!(b,2)
	facet!(b,p2,p3)
	facetregion!(b,3)
	facet!(b,p3,p4)	
	facetregion!(b,4)
	facet!(b,p4,p1)	
	
	simplexgrid(b,maxvolume=4.0^(-nref-1))
end


# ╔═╡ 84ad1220-2442-49c2-b33c-c5fe9e8d8b90
inpluto && gridplot(rect_t(nref=3,w=1.2))

# ╔═╡ d03965c5-8515-4fea-9a70-fe6431cfd6a4
function rect_r(;nref=0,w=1.0,h=1.0)
	X=0:2.0^(-nref-1):w
	Y=0:2.0^(-nref-1):h
	@assert X[end]≈w
	@assert Y[end]≈h
	simplexgrid(X,Y)
end

# ╔═╡ 61892958-943a-4531-8d28-278f6df011ff
inpluto && gridplot(rect_r(nref=3,h=1.25))

# ╔═╡ 79ba0ba5-2e34-44cb-beff-223d3e42cc2e
function hminmax(grid)
	en=grid[EdgeNodes]
    coord=grid[Coordinates]
	nedges=size(en,2)
	hmin=floatmax()
	hmax=-floatmax()
	for iedge=1:nedges
		@views d=norm(coord[:,en[1,iedge]]-coord[:,en[2,iedge]])
		hmin=min(hmin,d)
		hmax=max(hmax,d)
	end
	hmin,hmax
end

# ╔═╡ Cell order:
# ╟─582b68c7-bf18-47eb-a083-f58d2bdd60ab
# ╠═ba562efe-1157-11ec-3ed6-97e95cd04353
# ╟─9ca39b3c-fca9-4a53-adea-61caa83c05fd
# ╠═8e564f00-735d-4fca-821d-50fce8e84d79
# ╟─e2ca7b02-4233-4e72-ab26-0f0dcb7bd601
# ╠═6c0000d3-1ef8-43b4-86d6-e65a63577ae3
# ╠═408a3da2-0751-4e78-a8a3-d235f76ea13f
# ╟─87fa7938-4ee7-4749-9538-1116387b808f
# ╟─3b62a531-b4a2-4460-926d-03d869649139
# ╠═73a8bedc-09ee-4d81-adb4-6c02ca1ad86b
# ╟─32099f08-5181-4f60-b9cc-88fed3201a33
# ╠═8487f01e-6114-40c3-9d8a-92b3676d4c4c
# ╠═25c049e7-c9c4-4326-a501-34a16d48df7c
# ╠═626bbd8b-4818-4ad7-9924-0f274ba55561
# ╠═0deab116-197a-4cfc-a94a-529adfa4ce4f
# ╟─3375c6fe-0ffc-4a80-9f73-245d3f8a4abd
# ╠═92b4efbe-0bac-4a47-89e5-e482b3ea96be
# ╟─2a6e529c-cf62-482b-bc99-f0d705cf0ebd
# ╠═3db6ebb1-dc8c-4538-bf23-2c7d0b06ad09
# ╠═2878ae3f-8e4b-4f76-a194-d01b8993f6a8
# ╟─b9d2b069-fbd5-4d0a-a4b7-c87cdff5a03d
# ╠═b1c33203-6a8d-4748-87ba-9dfc237f90fc
# ╠═f58d0c1d-6009-4853-b288-5f1de2d90f40
# ╠═5bd9792b-04eb-4116-a85d-d6201acba562
# ╠═84ad1220-2442-49c2-b33c-c5fe9e8d8b90
# ╠═d03965c5-8515-4fea-9a70-fe6431cfd6a4
# ╠═61892958-943a-4531-8d28-278f6df011ff
# ╠═79ba0ba5-2e34-44cb-beff-223d3e42cc2e
