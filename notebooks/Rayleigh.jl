### A Pluto.jl notebook ###
# v0.19.21

using Markdown
using InteractiveUtils

# ╔═╡ 6482eb8c-453e-47d5-8488-95cb440bd3c1
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise
	using ExtendableSparse
	using SparseArrays
	using Printf
	using LinearAlgebra
	using Arpack
end

# ╔═╡ 9cc0927e-c32d-4ba4-85ea-70f975759e0d
md" ### Rayleigh quotient iteration
1. Initialize vector ``u_0``, ε tolerance \
2. ``\displaystyle\lambda_0 := \frac{u_0 ^* A u_0}{u_0^* M u_0}`` \
3. for ``k=1,2,..`` until convergence: ``\quad\displaystyle u_k := (A-\lambda_{k-1} MI)^{-1} u_{k-1}, \quad \lambda_k := \frac{u_k^* A u_k}{u_k^* M u_k}``"
		

# ╔═╡ f619f4d0-983e-401c-b773-f22004000255
function rayleigh_general(A,M,u,ε=1.0e-17)
	#check norms
	k = 0
	
	y=u
	ynorm = norm(y)   #oder ist hier die Norm sqrt(y'*M*y)?? Falscher eVec.
	#ynorm = sqrt(y'*M*y)
	y = y/ynorm
	λ = (u'*A*u) / (u'*M*u)
	
	u = (A - λ*M)\y
	unorm = norm(u)
	#unorm = sqrt(u'*M*u)
	u = u/unorm
	err = norm( u/unorm - (u'*y *y) / (ynorm^2 * unorm) )
	info=[]
	while err>ε && k<20
		k = k+1
		
		y = u
		ynorm = norm(y)
		#ynorm = sqrt(y'*M*y)
		#y = y/ynorm
		λ = (u'*A*u) / (u'*M*u)
		
		u = (A - λ*M)\y
		unorm = norm(u)
		#unorm = sqrt(u'*M*u)
		u = u/unorm
		
		err = norm( u/unorm - (u'*y *y) / (ynorm^2 * unorm) )
		
		push!(info,(k, λ, unorm, ynorm, err))
	end
	@info info	
	return (u,λ)
end

#überprüfe Voraussetzungen an den init Vektor


# ╔═╡ a45c84ca-a9d0-4755-aff5-eee9397a465b
function rayleigh_angle_crit(A, u, ε=1.0e-15) #for M=I
	#Abbruchkrit: Winkel zwischen approx. Vektoren (atkuelle und vorherige Iteration)
	k = 0
	unorm = norm(u)
	
	y = u/norm(u)
	λ = y' * A * y 
	
	u = (A - λ*I)\y
	unorm = norm(u)
	u = u / unorm
	
	err = norm(u - u'*y * y) #sin of angle between vectors
	
	info=[]
	while err > ε
		k = k+1
		y = u
		
		λ = u' * A * u 
		u = (A - λ*I)\u
		u = u/norm(u)
		
		err = norm(u - u'*y * y) #def.= y* y'*u, aber in Norm gleich
		push!(info, (k, λ, unorm, err))
	end
	@info info
	return (u,λ)
end

# ╔═╡ b1af7c63-5df3-4fa8-95f6-eff55746d2d1
function rayleigh(A,u) #for M=I
	#Abbruchkriterium: blow up der Norm des approx. Eigenvektors
	
	k = 0
	unorm = norm(u)
	
	u = u/norm(u)
	λ = u' * A * u 
	λ = 1.05 * λ
	info=[]
	while unorm <1.0e15 && k<20
		y = u
		ρ = λ
		
		k = k+1
		
		u = (A - λ*I)\u
		unorm = norm(u)
		if isfinite(unorm) == false
			return (y, ρ)
		else
		u = u/unorm
				
		push!(info,(k, λ, unorm))
		λ = u' * A * u 
		end
	end
	@info info	
	return (u,λ)
end

# ╔═╡ 173c6d0a-4e62-4e2c-83b9-3fa9a56d3e1b
function rayleigh_step(A, u) #for M=I
	k = 0
	unorm = norm(u)
	u = u/norm(u)
	λ = u' * A * u
	
		y = u
		ρ = λ
		
		
		λ = u' * A * u
		@info (A - λ*I), λ, u
	
		u = (A - λ*I)\u
		unorm = norm(u)
		u = u/unorm
		@info k, λ, unorm

		
	return (u,λ)
end

# ╔═╡ 7cfe6d84-453e-4b5d-8578-1db522977cce
A =[2 -1 0 0 0 0 0 0 0; -1 2 -1 0 0 0 0 0 0; 0 -1 2 -1 0 0 0 0 0; 0 0 -1 2 -1 0 0 0 0; 0 0 0 -1 2 -1 0 0 0 ;0 0 0 0 -1 2 -1 0 0; 0 0 0 0 0 -1 2 -1 0 ; 0 0 0 0 0 0 -1 2 -1; 0 0 0 0 0 0 0 -1 2] 

# ╔═╡ 81332cec-3919-4993-95cf-51b11afd8ac8
u = [-4, -3, -2, -1, 0, 1, 2, 3, 4]

# ╔═╡ 1b45bb12-c535-47c3-bf8a-2b286b6c9510
(vec2, lam2) = rayleigh_angle_crit(A,u)

# ╔═╡ 3b29dd01-39ba-4fae-8653-78e1b05ac95d
(vec1, lam1) = rayleigh(A,u) 

# ╔═╡ b19b3852-0bac-4b92-81d1-21b46d208ac1
(vek1, lamb1) = rayleigh_general(A,I,u)

# ╔═╡ 33d95f8f-0734-4b34-a56e-391ab2ab6903
D = Diagonal(1:5)

# ╔═╡ d6f36a4c-53b8-4d11-9df5-157f27570396
v = [1,2,3,0,-5]

# ╔═╡ 2de313b5-d52c-429d-a8f1-29ac8f3c3db3
H = [1 1; 1 3]

# ╔═╡ 1f6c5c91-1ffc-418d-905d-987162ea4ffe
z = [0,5]

# ╔═╡ ec94ed57-1bd9-4ca7-8f8c-f51c023ad350
(z1,λ)=rayleigh_step(H,[1,2])

# ╔═╡ af7bdd48-ec73-471f-9524-14e788e6ec14
(z2,λ2)=rayleigh_step(H, z1)

# ╔═╡ 97951a32-cbe1-48d7-b7c7-ff26ba6d94b9
(z3,λ3)=rayleigh_step(H, z2)

# ╔═╡ 7152ee9f-8251-4fb8-a893-9431808528c3
(z4,λ4)=rayleigh_step(H, z3)

# ╔═╡ 7d1ee348-6afc-4884-8767-8c95f7b21975
(z5,λ5)=rayleigh_step(H, z4)

# ╔═╡ 1b46be7a-e620-4b9e-94f2-5019917d1565
A1 = Diagonal(2:6)
	

# ╔═╡ 5ff6ceaf-8be6-45cb-9d24-e845bc0937a5
B1 = Diagonal(1:5)

# ╔═╡ 8061afe5-a577-4f4d-aa2d-17c9e0ad18ab
A2 = [1 2 3 4 0; 2 1 3 2 0; 3 3 2 0 4; 4 2 0 1 2; 0 0 4 2 1]

# ╔═╡ 749d25b3-0570-4285-8830-4e2d4ebb3fef
ρ, φ = eigs(A1,B1, nev=3)

# ╔═╡ db4789c5-11a7-4e3b-acc7-8774c0eee3bf
u1 = rand(5)

# ╔═╡ 2e1a36c2-b277-40e8-a27e-876520e211e7
rayleigh_general(A1,B1,u1)

# ╔═╡ 0b5384b1-5b67-42ba-8db4-d186a43f4a7a
ρ1, φ1 = eigs(A1,B1, nev=3)
#E.val for A1,B1:   2,   1.5,   1.3,   1.25,   1.2

# ╔═╡ 73424beb-3ad4-46a6-b991-587d37852eb1
rayleigh_general(A2,B1,u1)

# ╔═╡ 312597dd-d34d-4e1e-94d6-48c765e35fe7
ρ2, φ2 = eigs(A2,B1, nev=3)

# ╔═╡ Cell order:
# ╠═6482eb8c-453e-47d5-8488-95cb440bd3c1
# ╟─9cc0927e-c32d-4ba4-85ea-70f975759e0d
# ╠═f619f4d0-983e-401c-b773-f22004000255
# ╠═a45c84ca-a9d0-4755-aff5-eee9397a465b
# ╠═b1af7c63-5df3-4fa8-95f6-eff55746d2d1
# ╠═173c6d0a-4e62-4e2c-83b9-3fa9a56d3e1b
# ╠═1b45bb12-c535-47c3-bf8a-2b286b6c9510
# ╠═3b29dd01-39ba-4fae-8653-78e1b05ac95d
# ╠═b19b3852-0bac-4b92-81d1-21b46d208ac1
# ╟─7cfe6d84-453e-4b5d-8578-1db522977cce
# ╟─81332cec-3919-4993-95cf-51b11afd8ac8
# ╟─33d95f8f-0734-4b34-a56e-391ab2ab6903
# ╟─d6f36a4c-53b8-4d11-9df5-157f27570396
# ╟─2de313b5-d52c-429d-a8f1-29ac8f3c3db3
# ╟─1f6c5c91-1ffc-418d-905d-987162ea4ffe
# ╠═ec94ed57-1bd9-4ca7-8f8c-f51c023ad350
# ╠═af7bdd48-ec73-471f-9524-14e788e6ec14
# ╠═97951a32-cbe1-48d7-b7c7-ff26ba6d94b9
# ╠═7152ee9f-8251-4fb8-a893-9431808528c3
# ╠═7d1ee348-6afc-4884-8767-8c95f7b21975
# ╟─1b46be7a-e620-4b9e-94f2-5019917d1565
# ╟─5ff6ceaf-8be6-45cb-9d24-e845bc0937a5
# ╟─8061afe5-a577-4f4d-aa2d-17c9e0ad18ab
# ╠═749d25b3-0570-4285-8830-4e2d4ebb3fef
# ╠═db4789c5-11a7-4e3b-acc7-8774c0eee3bf
# ╠═2e1a36c2-b277-40e8-a27e-876520e211e7
# ╠═0b5384b1-5b67-42ba-8db4-d186a43f4a7a
# ╠═73424beb-3ad4-46a6-b991-587d37852eb1
# ╠═312597dd-d34d-4e1e-94d6-48c765e35fe7
