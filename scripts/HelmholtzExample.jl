module HelmholtzExample

using ExtendableGrids
using GridVisualize
using SimplexGridFactory
using Triangulate
using ExtendableSparse
using HelmholtzProject




function test1d(;N=500,kwargs...)
    X=0:1/N:1
    grid=simplexgrid(X)
    λ,v=eigensolve(grid;kwargs...)
    grid,λ,v
end

function test2d(;N=50,kwargs...)
    X=0:1/N:1
    grid=simplexgrid(X,X)
    λ,v=eigensolve(grid;kwargs...)
    grid,λ,v
end

function test1d_confined(;N=500,k=200+1im,kwargs...)
    X=0:1/N:1
    grid=simplexgrid(X)
    λ,v=eigensolve(grid;
                   n=(i,x)->( abs(x[1]-1/2) < 0.1 ? 1 : 0.5 ),
                   k=k,
                   kwargs...)
    grid,λ,v
end


function test2d_confined(;nref=4,k=200+1im,kwargs...)
    grid=HelmholtzProject.grid2d_circle_in_square(nref=nref)
    n(i,coord)= i==1 ? 0.5 : 1 
    λ,v=eigensolve(grid;n,k,kwargs...)
    grid,λ,v
end


#=
function run1d(;N=500, n=(x)->( abs(x-1/2) < 0.1 ? 1 : 0.5 ))
    X=0:1/N:1
    grid=simplexgrid(X)
    



# ╔═╡ b9733cd3-502d-4ff1-912b-a02bfd1fd333
b1(i,coord)= (0, 0) 

# ╔═╡ 1ffd1959-97ee-4c18-8607-93917b28300b
md"""
real(k)= $(@bind rkT Slider(0:0.5:500,show_value=true,default=200))
imag(k)= $(@bind ikT Slider(-50:0.1:50,show_value=true,default=1.0))
"""

# ╔═╡ d1a2d08e-2519-4c2e-b9a8-6a07fd2ad844
kT=rkT+ikT*im;kT^2

# ╔═╡ 6727d195-5897-40d9-9a66-a6e572736783
viseT1=GridVisualizer(resolution=(600,150),legend=:rt); viseT1

# ╔═╡ 16d201b6-b352-4130-a9f7-4dd702319442
viseT2=GridVisualizer(resolution=(600,150),legend=:rt); viseT2
λ=
# ╔═╡ 2ab796c4-ad9c-469f-a358-6aa48fbcfba0
function plotev(λ)
	vis=GridVisualizer(size=(500,200),legend=:rt)
	r=1:length(λ)
	scalarplot!(vis,r,abs.(λ),label="|λ|",color=:green)
	scalarplot!(vis,r,real.(λ),label="Re(λ)",color=:red,clear=false)
	scalarplot!(vis,r,imag.(λ),label="Im(λ)",color=:blue,clear=false)
	reveal(vis)
end

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
@time λ1,v1=jd_eigen(A1,M1,nev=20);

# ╔═╡ ddca0556-9b43-4125-989d-e1a0547a3806
md"""
Eigenvalue index $(@bind iT PlutoUI.Slider(1:length(λ1),show_value=true))
"""

# ╔═╡ 52111b28-fc99-413c-a702-88bc252f383e
md"""λ($iT) = $(round(λ1[iT],sigdigits=5))"""

# ╔═╡ 089c2f36-fe58-4450-8aed-cb3f7825b1a2
begin
	z=0.0-1.0im
    @views uet=v1[:,iT]*z/abs(z)
    s=sign(real(uet[2]))
    scalarplot!(viseT1,
        X1,s*imag(uet),label="im",color=:blue)
    scalarplot!(viseT1,
        X1,real(uet)*s,label="re",clear=false,color=:red,show=true)
end

# ╔═╡ e3edbbd9-19e4-41ee-8cba-05575b05c1b1
scalarplot!(viseT2,
        X1,abs.(uet),label="abs",color=:green,show=true)

# ╔═╡ 5d18de00-a949-49a7-9a9e-650da7c2bdbb
plotev(λ1)

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
@time λ,v=arpack_eigen(A2,M2;nev=nev,maxiter=500);

# ╔═╡ ba542ed8-5408-4b47-9877-0d6651aa472c
md"""
Eigenvalue index $(@bind iT2 PlutoUI.Slider(1:length(λ),show_value=true))
"""

# ╔═╡ 5aec611b-94a9-420b-a71c-64b09ee7e373
λ[iT2]

# ╔═╡ 9dee97de-c7c8-4209-a2bd-c22a3e394b67
vis=GridVisualizer(dim=2);vis

# ╔═╡ db8521f3-9ceb-4f63-8032-9da454d54577
scalarplot!(vis,grid2,abs.(v[:,iT2]),label="abs",colormap=:hot,show=true,levels=0)

# ╔═╡ 2e873ab7-08e1-4d8b-acbe-dcaec4218ba9
plotev(λ)
    =#
end
