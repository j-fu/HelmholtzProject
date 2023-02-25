### A Pluto.jl notebook ###
# v0.19.22

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

# ╔═╡ 38fc344e-ebe2-11eb-0486-af4c3d700950
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using VoronoiFVMDiffEq
	using DifferentialEquations
	using ExtendableGrids
	using GridVisualize
	using PlutoUI
	using PlutoVista
	using ExtendableSparse
	using LinearAlgebra
	using Logging
	using TerminalLoggers
	default_plotter!(PlutoVista)
end

# ╔═╡ 82ce69ac-75d4-410a-a237-ed466af0c762
md"""
# Wave and Helmholtz equations
"""

# ╔═╡ a65b66e2-3cf4-4df5-91ce-fd45a920899d
TableOfContents(depth=4,aside=true)

# ╔═╡ a5a753a5-450b-42d1-9511-7d39951bf9db
md"""
## The wave equation as system of equations
"""

# ╔═╡ b5287fbe-c149-4b4c-87ae-c095188fb06e
md"""
This is the n-dimensional wave equation:

```math
u_{tt}- c^2 \Delta u = 0
```


We can create a system of first oder in time PDEs out of this:

```math 
    \begin{aligned}
        u_t - v&=0\\
     v_t -c^2\Delta u&=0
\end{aligned}
```

This allows for a quick implementation in VoronoiFVM (which may be not the optimal way, in particular with respect to time discretization).
"""

# ╔═╡ b354f875-9f6c-418e-add3-c34f85f1d23e
const iu=1; const iv=2;

# ╔═╡ 37fa4962-6b74-4c92-b716-8c1d2e9cfc00
storage(y,u,node,data)=y.=u;

# ╔═╡ a9febfb8-05b4-4bd6-828d-ddccc3c9bdc9
reaction(y,u,node,data)= y[iu]=-u[iv];

# ╔═╡ 1e01a428-4dc5-4a57-877f-0c9d84e70149
flux(y,u,edge,data)=y[iv]=data.c^2*(u[iu,1]-u[iu,2]);

# ╔═╡ a4241f4b-6ce9-4024-8145-5d94fe8fcbbc
md"""
Wave velocity: 
"""

# ╔═╡ 6ebb83f7-4b98-40fd-be1b-ebb2ce4399a0
const c=0.1

# ╔═╡ 779fbc09-5801-4c2d-8a74-e0ea576af121
md"""
Domain length:
"""

# ╔═╡ f4f2a727-8463-4845-b085-55a42be8da51
L=4

# ╔═╡ fac514b1-8088-4031-9f61-6851aeb916bc
N=151

# ╔═╡ ccdcf06e-3026-43e4-9624-65f36c232605
md"""
### Perturbation in the center of the domain 
"""

# ╔═╡ 5c903d43-1caa-473d-9fc2-f11c13b70cb8
grid=simplexgrid(range(-L,L,length=N));

# ╔═╡ 196ac156-9999-4e0f-92c0-2387d84b3555
dt=1.0e-2; tend=100.0;

# ╔═╡ 1781ac71-de97-48c9-94fa-a4b59c74a7e3
md"""
Implementation of transparent or mirror bc
"""

# ╔═╡ 707c0d50-065a-44fd-a3a6-33210b8bcf90
function brea(y,u,node,data)
	if node.region==2 
		if data.bctype==:transparent
	     	y[iu]=data.c*u[iu]
		elseif data.bctype==:mirror
			boundary_dirichlet!(y,u,node,species=1,region=2,value=0)
		end
	end
end;

# ╔═╡ 50650c99-3ed8-48b0-86b0-8e2d63bfb797
md"""
Package wave number κ: $(@bind κ Scrubbable(0:5,default=0))
"""

# ╔═╡ c61ad8fe-a8b8-4dd8-8852-45a41c5884b0
md"""
Boundary condition at x=L: $(@bind bc2type Select(["reflection","mirror","transparent"]))
"""

# ╔═╡ 7d3c2f93-400d-426b-b00d-c7d84ab02633
sys=VoronoiFVM.System(grid,storage=storage,flux=flux,breaction=brea, reaction=reaction,data=(c=c,bctype=Symbol(bc2type)), species=[1,2])


# ╔═╡ 8fb23cce-7c26-4b62-b03e-39044551caa8
begin
   inival=unknowns(sys,inival=0)	
   inival[1,:].=map(x->cos(κ*π*x)*exp(-x^2/0.1) ,grid)	
end;

# ╔═╡ 00b362c3-eaec-4e9e-a119-bf22cd7aa33e
problem = ODEProblem(sys,inival,(0.0,tend))

# ╔═╡ 4f14fbc5-e369-497d-b234-37460809a0bd
tsol=DifferentialEquations.solve(problem,
	                               Rosenbrock23();  
		                           initializealg=NoInit(),
                                   force_dtmin=true,
                                   adaptive=true,
                                   reltol=1.0e-4,
                                   abstol=1.0e-5,
                                   dtmin=dt,
	                               progress=true,
	                               progress_steps=1,
                                   dt=dt);

# ╔═╡ 68fefe0a-caf6-4ace-9119-887e12f85c17
tsol.alg

# ╔═╡ 6772197f-a8fb-44bb-8ffb-294dba34c1fc
tsol1=reshape(tsol,sys);

# ╔═╡ ff86aa0c-108a-46fb-927d-f4ce4423e17c
if bc2type=="reflection"
	md"""Reflection (Neumann) bc ``\partial_x u|_{x=L}=0``"""
elseif bc2type=="mirror"
	md"""Mirror (Dirichlet) bc ``u|_{x=L}=0``"""
elseif bc2type=="transparent"
	md"""Transparent  bc ``(\partial_t u - c\partial_xu)|_{x=0}=0`` """
end


# ╔═╡ 98b03f91-b0ac-4fa2-b0f3-62075500b43a
md"""
t=$(@bind t PlutoUI.Slider(range(0,tsol.t[end],length=10001);show_value=true))
"""

# ╔═╡ d09bb020-0504-465d-91bd-f7e7453ced2a
vis=GridVisualizer(resolution=(600,150));vis

# ╔═╡ 63f27b01-4f9a-44db-95ce-d636e67a309c
let
	u=tsol1(t)
	scalarplot!(vis,grid,u[1,:],flimits=(-1,1),clear=true,show=true,title="t=$(t)")
	vis
end

# ╔═╡ 0c2fd4f9-8164-46df-9990-b74440e01c2a
md"""
### Outer excitation
"""

# ╔═╡ 52a5f159-4474-49d0-905f-dda0abc25a43
md"""
Here, we set time periodic Dirichlet bc at x=0: ``u|_{x=0}=\sin(ω t)`` and transparent BC at x=L
"""

# ╔═╡ da5d6b61-f1e4-4aec-85cb-91503cbf2732
function brea2(y,u,node,data)
	if node.region==1
		y[iu]=1.0e5*(u[iu]-sin(data.ω*node.time))
	elseif node.region==2 
	    y[iu]=data.c*u[iu]
	end
end;

# ╔═╡ 1c3e030a-a7eb-4854-aca6-a23d4cb783e8
ω=π/4

# ╔═╡ 1c93b9cb-8df3-4771-b828-061406100840
begin
	grid2=simplexgrid(range(0,2*π,length=301))
	sys2=VoronoiFVM.System(grid2,storage=storage,
		                         flux=flux,
		                         reaction=reaction,
		                         data=(c=c,ω=ω),
		                         breaction=brea2,
		                         species=[1,2])
   inival2=unknowns(sys2,inival=0)
end

# ╔═╡ 6d01f814-5944-49e0-a691-ac43061e91c3
begin
	dt2=1.0e-2
	tend2=100
	control2=VoronoiFVM.NewtonControl()
    control2.force_first_step=true
	control2.Δt_min=dt2
	control2.Δt_max=dt2
	control2.Δu_opt=100
	control2.verbose=false
end;

# ╔═╡ 478776ea-df69-4378-bf50-ab3f56335ffc
begin
	problem2 = ODEProblem(sys2,inival2,(0.01,tend2))
 	tsol20=DifferentialEquations.solve(problem2,
	                               Rosenbrock23(); 
		                           initializealg=NoInit(),
                                   force_dtmin=false,
                                   adaptive=true,
                                   reltol=1.0e-2,
                                   abstol=1.0e-3,
                                   dtmin=0.001*dt2,
	                               progress=true,
	                               progress_steps=1,
                                   dt=dt2);
	
tsol2=reshape(tsol20,sys2);
end

# ╔═╡ 9fd39ff9-e2c7-4dd0-831f-2f904706d2c5
md"""
t=$(@bind t2 PlutoUI.Slider(range(tsol2.t[1],tsol2.t[end],length=1001)))
"""

# ╔═╡ b10270aa-e90a-43ea-883b-79eb0b9ae105
vis2=GridVisualizer(resolution=(600,150),xlabel="x",ylabel="y",legend=:lt);vis2

# ╔═╡ a132a086-cd12-49ed-95e3-bf0edb71af9d
let
	u=tsol2(t2)[1,:]
	scalarplot!(vis2,grid2,u,flimits=(-2,2),show=true,clear=true,label="t=$(t2)")
end

# ╔═╡ a9fddbe6-7a44-4e75-a14d-62cfb9eb2ecf
md"""
## Helmholtz equation
"""

# ╔═╡ 460a2d8a-0619-417e-b442-30b326743ac1
md"""
### Longitudinal Helmholtz equation
"""

# ╔═╡ 0b0f1390-f4d1-4adb-be44-06dcf3b1c26f
md"""

#### Derivation by separation of variables
Let us write the wave equation as

```math
u_{zz}= \frac{1}{c^2 }u_{tt}
```

and assume that ``u(z,t)= ζ(z)\psi(t)``

Then the equation transforms to
```math
\zeta_{zz}\psi= \frac{1}{c^2 }\zeta \psi_{tt}
```

and, by assuming a constant ``k``

```math
\frac{\zeta_{zz}}{\zeta}=-k^2=\frac{1}{c^2 }\frac{\psi_{tt}}{\psi}
```

This gives 
```math
\begin{aligned}
  \zeta_{zz} + k^2 \zeta &=0\\
  \psi_{tt} + ω^2 \psi &=0
\end{aligned}
```

with ``ω=kc``.



"""

# ╔═╡ 66eb2d95-169a-4933-9d96-3589a680c05a
md"""
#### Derivation by  time periodic separation ansatz
```math
\begin{aligned}
  u(z,t)&=\zeta(x)e^{i\omega t}\\
  u_{zz}&=\zeta_{zz}e^{i\omega t}\\
  u_{zz}&=-\omega^2\zeta e^{i\omega t}\\
  \zeta_{zz}e^{i\omega t}&=-\frac{\omega^2}{c^2}\zeta e^{i\omega t}\\
  \zeta_{zz}+k^2\zeta&=0
\end{aligned}
```

with ``k=\frac{\omega}{c}``

Transparent boundary conditions translate as follows:

```math
\begin{aligned}
\partial_t u &=iω\zeta e^{i\omega t}\\
\partial_z u &= \zeta_z e^{i\omega t}\\
\partial_t u - c\partial_z u&=0\\
i\omega	 ζ - c\zeta_z&=0\\
-\zeta_z + ik  ζ&=0
\end{aligned}
```
"""

# ╔═╡ 24d8f64a-5e5b-4861-9744-2450ea557c60
md"""
#### Numerical solution

Once again, time periodic Dirichlet bc at x=0 and reflection BC at x=L.
Time periodic Dirichlet bc translate into  ``\phi|_{x=0}=1``
"""

# ╔═╡ d5fb9a17-f37b-4138-8a2d-55756013d639
md"""
Excitation at "inlet" translates as ``\zeta|_{z=0}=1``
"""

# ╔═╡ a8b323c3-ee98-438c-85cb-dc81f4445f73
N1=201

# ╔═╡ fdbb9b0b-164e-4e0d-89f8-759ca9fd9fe7
X=collect(range(0,1,length=N1));

# ╔═╡ 6e141435-4875-4747-ab8c-b24d18f0f423
ζ0=1;

# ╔═╡ e00491e2-41cb-4a1b-8b1b-b5fed67c1802
function assemble(T,N1,k;penalty=false,transparent=false)
    A=ExtendableSparseMatrix(T,N1,N1)
    b=zeros(T,N1) # rhs
    m=zeros(T,N1) # -> mass matrix
    for i=1:N1-1 # loop over all intervals (X[i],X[i+1])."Voronoi cells" are 
    # Ω_i=( (X[i-1]+X[i])/2, (X[i]+X[i+1])/2)
        j=i+1
        h=X[j]-X[i]
        # "stiffness matrix" - ζ_{zz}    
        A[i,i]+=1/h
        A[i,j]+=-1/h
        A[j,j]+=1/h
        A[j,i]-=1/h
        
        A[i,i]-=k^2*h/2 # Anteil an Ω_i    
        A[j,j]-=k^2*h/2 # Anteil an Ω_{i+1}    
        
        m[i]+=h/2
    m[j]+=h/2    
    end
    # Transparent bc at z=L 
    if transparent
        A[end,end]+=k*1im
    end
    
    if penalty
        pen=1.0e10
        A[1,1]+=pen
        b[1]=ζ0*pen
        m[1]=1
    else # this gives spurios eigensolutions
        A[1,1]=1
        A[1,2]=0
        b[1]=1
        m[1]=1        
    end
    
    M=Diagonal(m)
    A,M,b
end;

# ╔═╡ 97b3ca06-ead1-4a08-b955-d3604cba1250
md"""
Solve `` A(k)ζ=b`` where ``A(k)=-Δ_h-k^2M`` 
"""

# ╔═╡ ae37271b-ea3f-418d-8a1e-36292bc65b6b
md"""
Penalty: $(@bind penalty CheckBox())
Transparent: $(@bind transparent CheckBox())
"""

# ╔═╡ 0ebd628e-7561-4a1a-8374-bf61d8b6675c
md"""
real(k)= $(@bind rk Slider(0:0.1:100,show_value=true))
imag(k)= $(@bind ik Slider(-50:0.1:50,show_value=true,default=0.0))
"""

# ╔═╡ 7d1338f9-6f63-4083-bcce-dac6c785682a
k=rk+ik*im; md""" k^2= $(k^2)"""

# ╔═╡ 1701fe04-640d-43e2-807a-9fd7ae9e5ad5
A,M,b=assemble(Complex{Float64},N1,k,penalty=penalty,transparent=transparent);

# ╔═╡ 605bbd80-d710-4499-b1b8-038c0bb478ba
ζ=A\b;

# ╔═╡ a62585d0-41f4-495c-87ee-7d96d78b5d61
vis3=GridVisualizer(resolution=(600,150),legend=:rt); vis3

# ╔═╡ a240a6a1-00bf-4239-b099-d5f8d55322b6
vis4=GridVisualizer(resolution=(600,150),legend=:rt); vis4

# ╔═╡ 13d8eec4-1cd9-4c90-8228-7ee8b6b667bc
begin
scalarplot!(vis3,X,real(ζ),label="re ζ",color=:red,clear=true)
scalarplot!(vis3,X,imag(ζ),label="im ζ",color=:blue,clear=false,show=true)
end

# ╔═╡ 23a0afb0-f571-406f-98ab-0bca898a83fe
visev=GridVisualizer(resolution=(600,150),legend=:rt); visev

# ╔═╡ db238b44-9cd7-4cef-9391-c9040687e52c
scalarplot!(vis4,X,abs.(ζ),label="|ζ|",color=:green,flimits=(0,2),show=true)

# ╔═╡ 2bb3d025-bf3f-4c1c-9b25-3f7032bbd2a2
md"""
Compute eigenvalues ``-\Delta_h \zeta - k^2M\zeta = \lambda M\zeta``: $(@bind compute_ev PlutoUI.CheckBox())
Eigenvalue index $(@bind idx PlutoUI.Slider(1:N1,show_value=true))
"""

# ╔═╡ 9dc672df-4c52-476c-8863-a4bdc716836d
if compute_ev
    ev=eigen(Matrix(A),Matrix(M))
end;

# ╔═╡ 4c4a28fe-be0d-4391-8619-00b7ea9760a0
if compute_ev
md"""
λ($(idx)) = $(round(ev.values[idx],sigdigits=5))
"""
end

# ╔═╡ fb635d85-7751-4126-aad4-4eb1558a0098
if compute_ev let
    u=ev.vectors[:,idx]
    s=sign(real(u[2]))
    remax=maximum(abs.(real(u)))
    imax=maximum(abs.(imag(u)))
        title=""

    scalarplot!(visev,X,s*imag(u),label="im",color=:blue,title=title,flimits=(-2,2))
    scalarplot!(visev,X,real(u)*s,label="re",color=:red,clear=false,show=true,title=title,flimits=(-2,2))
end
end

# ╔═╡ ec866adb-8764-4cc9-90fb-51af7ad64b1b
md"""
- We see a significant increase in ``|\zeta|`` if ``k^2`` is close to an eigenvalue for the Problem with ``k=0`` due to resonance
- The introductionn of Dirichlet bounday conditions in the way performed here introduces spurious eigenfunctions and eigenvalues:
   - For the penalty method, we get a largest eigenvalue approximately of the value of the penalty parameter
   - For the matrix modification method, we get an eigenvalue 1 and an eigenfuntion which is nonzero at ``x=0``.  This value is in the order of magnitude of other eigenvalues and would be hard to "filter out". Note that we could scale the "Dirichlet equation" for `u[1]`, shifting this value to a range where it does not interfere with the relevant eigenvalues.
"""

# ╔═╡ 182cefbf-27c6-4c6e-93a4-6066b08fc356
md"""
### Transversal-longitudinal separation
"""

# ╔═╡ a5ac0e28-a93b-4ff1-b0b9-71a008f2ba9c
md"""
``\Delta u = \frac{1}{c^2} u_{tt}``

Ansatz:
```math
\begin{align}
  u(x,y,z,t) &= \phi(x,y)e^{i\kappa z}e^{i\omega t}\\
  \Delta u &= (\Delta_{xy}\phi  - k_0^2\phi)e^{i\kappa z}e^{i\omega t}\\
  \Delta_{xy}\phi + (\frac{\omega^2}{c^2}- \kappa^2)\phi &= 0
\end{align}
```
"""

# ╔═╡ d35fe4d2-1d06-4871-8b65-57678de3d484
md"""
### Numerical solution of 1D transversal problem

Here, we want to use transparent boundary conditions which will be not
possible only approximately (if at all) in 2D.

In addition, we get closer  to the laser problem by regarding the eigenvalue problem
``
\phi_{xx} + n(x)k^2\phi = \lambda\phi
``
which comes from a space dependent refractive index (wave speed over speed of light) ``n(x)``.
"""

# ╔═╡ 3db3df59-0609-40f9-a589-21c9a7269696
function assembleT(T,N1,n,k;transparent=false)
    A=ExtendableSparseMatrix(T,N1,N1)
    b=zeros(T,N1)
    m=zeros(T,N1)
    for i=1:N1-1
        j=i+1
        h=X[j]-X[i]
        A[i,i]+=1/h - n(X[i])*k^2*h/2
        A[i,j]+=-1/h
        A[j,j]+=1/h - n(X[j])*k^2*h/2
        A[j,i]-=1/h
        m[i]+=h/2
        m[j]+=h/2    
    end
    if transparent
        A[end,end]+=k*1im
        A[1,1]+=k*1im
    end
    
    M=Diagonal(m)
    A,M,b
end;

# ╔═╡ c9c1d9e4-d793-4cd9-8dbe-128831f6bb17
n(x)= abs(x-1/2) < 0.1 ? 1 : 0.5 ;

# ╔═╡ 04229b30-3dc9-4c34-b1af-c69e095f3d2b
scalarplot(X,n.(X),resolution=(600,150))

# ╔═╡ 96a96a52-919d-4191-9c4c-93624445dbc5
md"""
Transparent: $(@bind transparentT CheckBox())
real(k)= $(@bind rkT Slider(0:0.5:500,show_value=true))
imag(k)= $(@bind ikT Slider(-50:0.1:50,show_value=true,default=0.0))
"""

# ╔═╡ d80072d1-b9d2-4141-acb0-61d9040456e9
kT=rkT+ikT*im;kT^2

# ╔═╡ 036dce1e-8537-4ad1-835d-fdd89cdcfdbf
AT,MT,bT=assembleT(Complex{Float64},N1,n,kT,transparent=transparentT);

# ╔═╡ 1b49f27c-215f-4a4a-a690-50f8317873bc
eT=eigen(Matrix(AT),Matrix(MT));

# ╔═╡ 1086b569-6549-4453-9d33-61c77bc76b75
md"""
Eigenvalue index $(@bind iT PlutoUI.Slider(1:N1,show_value=true))
"""

# ╔═╡ 83f897a5-751e-45e5-8e8c-b0c8abf7e456
md"""λ($iT) = $(round(eT.values[iT],sigdigits=5))"""

# ╔═╡ 894c961b-42a6-4c86-be06-1e24e0306839
viseT1=GridVisualizer(resolution=(600,150),legend=:rt); viseT1

# ╔═╡ 659b11f3-7edf-43cc-a250-37a7344d82ff
viseT2=GridVisualizer(resolution=(600,150),legend=:rt); viseT2

# ╔═╡ c576898d-5304-4533-9493-b0075a878829
begin
    uet=eT.vectors[:,iT]
    s=sign(real(uet[2]))
    scalarplot!(viseT1,
        X,s*imag(uet),label="im",color=:blue,flimits=(-3,3))
    scalarplot!(viseT1,
        X,real(uet)*s,label="re",clear=false,color=:red,flimits=(-3,3),show=true)
end

# ╔═╡ 298ce775-d8b2-4d6f-b777-ada8ec9c4099
scalarplot!(viseT2,
        X,abs.(uet),label="abs",color=:green,flimits=(0,3),show=true)


# ╔═╡ 8c36312f-51b7-4490-b2b9-cd4561ba1f2e
md"""
- For a given eigenvalue index, we see for increasing ``k`` that the corresponding eigenfunction becomes confined in the region of large ``n``.
"""

# ╔═╡ Cell order:
# ╟─82ce69ac-75d4-410a-a237-ed466af0c762
# ╠═38fc344e-ebe2-11eb-0486-af4c3d700950
# ╠═a65b66e2-3cf4-4df5-91ce-fd45a920899d
# ╟─a5a753a5-450b-42d1-9511-7d39951bf9db
# ╟─b5287fbe-c149-4b4c-87ae-c095188fb06e
# ╠═b354f875-9f6c-418e-add3-c34f85f1d23e
# ╠═37fa4962-6b74-4c92-b716-8c1d2e9cfc00
# ╠═a9febfb8-05b4-4bd6-828d-ddccc3c9bdc9
# ╠═1e01a428-4dc5-4a57-877f-0c9d84e70149
# ╟─a4241f4b-6ce9-4024-8145-5d94fe8fcbbc
# ╠═6ebb83f7-4b98-40fd-be1b-ebb2ce4399a0
# ╟─779fbc09-5801-4c2d-8a74-e0ea576af121
# ╠═f4f2a727-8463-4845-b085-55a42be8da51
# ╠═fac514b1-8088-4031-9f61-6851aeb916bc
# ╟─ccdcf06e-3026-43e4-9624-65f36c232605
# ╠═5c903d43-1caa-473d-9fc2-f11c13b70cb8
# ╠═7d3c2f93-400d-426b-b00d-c7d84ab02633
# ╠═8fb23cce-7c26-4b62-b03e-39044551caa8
# ╠═196ac156-9999-4e0f-92c0-2387d84b3555
# ╠═00b362c3-eaec-4e9e-a119-bf22cd7aa33e
# ╠═4f14fbc5-e369-497d-b234-37460809a0bd
# ╠═68fefe0a-caf6-4ace-9119-887e12f85c17
# ╠═6772197f-a8fb-44bb-8ffb-294dba34c1fc
# ╟─1781ac71-de97-48c9-94fa-a4b59c74a7e3
# ╠═707c0d50-065a-44fd-a3a6-33210b8bcf90
# ╟─50650c99-3ed8-48b0-86b0-8e2d63bfb797
# ╟─c61ad8fe-a8b8-4dd8-8852-45a41c5884b0
# ╟─ff86aa0c-108a-46fb-927d-f4ce4423e17c
# ╟─98b03f91-b0ac-4fa2-b0f3-62075500b43a
# ╟─d09bb020-0504-465d-91bd-f7e7453ced2a
# ╠═63f27b01-4f9a-44db-95ce-d636e67a309c
# ╟─0c2fd4f9-8164-46df-9990-b74440e01c2a
# ╟─52a5f159-4474-49d0-905f-dda0abc25a43
# ╠═da5d6b61-f1e4-4aec-85cb-91503cbf2732
# ╠═1c93b9cb-8df3-4771-b828-061406100840
# ╠═1c3e030a-a7eb-4854-aca6-a23d4cb783e8
# ╠═6d01f814-5944-49e0-a691-ac43061e91c3
# ╠═478776ea-df69-4378-bf50-ab3f56335ffc
# ╟─9fd39ff9-e2c7-4dd0-831f-2f904706d2c5
# ╟─b10270aa-e90a-43ea-883b-79eb0b9ae105
# ╟─a132a086-cd12-49ed-95e3-bf0edb71af9d
# ╟─a9fddbe6-7a44-4e75-a14d-62cfb9eb2ecf
# ╟─460a2d8a-0619-417e-b442-30b326743ac1
# ╟─0b0f1390-f4d1-4adb-be44-06dcf3b1c26f
# ╟─66eb2d95-169a-4933-9d96-3589a680c05a
# ╟─24d8f64a-5e5b-4861-9744-2450ea557c60
# ╟─d5fb9a17-f37b-4138-8a2d-55756013d639
# ╠═fdbb9b0b-164e-4e0d-89f8-759ca9fd9fe7
# ╠═a8b323c3-ee98-438c-85cb-dc81f4445f73
# ╠═6e141435-4875-4747-ab8c-b24d18f0f423
# ╠═e00491e2-41cb-4a1b-8b1b-b5fed67c1802
# ╟─97b3ca06-ead1-4a08-b955-d3604cba1250
# ╠═1701fe04-640d-43e2-807a-9fd7ae9e5ad5
# ╠═605bbd80-d710-4499-b1b8-038c0bb478ba
# ╟─ae37271b-ea3f-418d-8a1e-36292bc65b6b
# ╟─7d1338f9-6f63-4083-bcce-dac6c785682a
# ╟─0ebd628e-7561-4a1a-8374-bf61d8b6675c
# ╟─a62585d0-41f4-495c-87ee-7d96d78b5d61
# ╟─a240a6a1-00bf-4239-b099-d5f8d55322b6
# ╟─13d8eec4-1cd9-4c90-8228-7ee8b6b667bc
# ╟─4c4a28fe-be0d-4391-8619-00b7ea9760a0
# ╟─23a0afb0-f571-406f-98ab-0bca898a83fe
# ╟─db238b44-9cd7-4cef-9391-c9040687e52c
# ╟─2bb3d025-bf3f-4c1c-9b25-3f7032bbd2a2
# ╟─fb635d85-7751-4126-aad4-4eb1558a0098
# ╟─9dc672df-4c52-476c-8863-a4bdc716836d
# ╟─ec866adb-8764-4cc9-90fb-51af7ad64b1b
# ╟─182cefbf-27c6-4c6e-93a4-6066b08fc356
# ╟─a5ac0e28-a93b-4ff1-b0b9-71a008f2ba9c
# ╟─d35fe4d2-1d06-4871-8b65-57678de3d484
# ╠═3db3df59-0609-40f9-a589-21c9a7269696
# ╠═c9c1d9e4-d793-4cd9-8dbe-128831f6bb17
# ╠═04229b30-3dc9-4c34-b1af-c69e095f3d2b
# ╠═036dce1e-8537-4ad1-835d-fdd89cdcfdbf
# ╠═d80072d1-b9d2-4141-acb0-61d9040456e9
# ╟─96a96a52-919d-4191-9c4c-93624445dbc5
# ╠═1b49f27c-215f-4a4a-a690-50f8317873bc
# ╟─1086b569-6549-4453-9d33-61c77bc76b75
# ╟─83f897a5-751e-45e5-8e8c-b0c8abf7e456
# ╠═894c961b-42a6-4c86-be06-1e24e0306839
# ╠═659b11f3-7edf-43cc-a250-37a7344d82ff
# ╠═c576898d-5304-4533-9493-b0075a878829
# ╠═298ce775-d8b2-4d6f-b777-ada8ec9c4099
# ╟─8c36312f-51b7-4490-b2b9-cd4561ba1f2e
