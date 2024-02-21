using DifferentialEquations, LinearAlgebra, LaTeXStrings, ProgressMeter, Statistics, Distributions
include("lattice.jl")

printstyled("VEGETATION TURBIDITY LATTICE MODEL\n"; bold=true, underline=true, color=:light_blue)

# Define the lattice
n = 20
m = 20
N = n*m
L = Lattice(n, m)

# Initial condition of the lattice
x0 = 1.0::Float64.*ones(Float64, length(L.grid))
# Initial condition of the bifurcation parameter (degree of smooth muscle activation)
c = 3.0::Float64 # background turbidity
push!(x0, c)
# SDE's parameters
σ = 0.05::Float64 # noise level
push!(x0, σ)
ε = 1.00::Float64 # parameter's drift
push!(x0, ε)

# Model's parameters
R = 0.2::Float64 # dispersion rate
rv = 0.5::Float64 # maximum growth rate 
hv = 0.2::Float64 # half-saturation cover
r = rand(Uniform(0.6::Float64, 1.0::Float64), N) # half-saturation turbidity

# Sliding window's parameters
T = 01.50
δt = 1e-3
width = 20

# Definition of the SDE - deterministic part
function iip_det!(f, x, neigh, t)
        for j=1:N
                E = hv*x[N+1]/(hv + x[j])
                T = (r[j]^4 + E^4)/(r[j]^4) 
                f[j] = rv*x[j]*(1.0::Float64-x[j]*T) + R*(x[neigh[1,j]] + x[neigh[2,j]] + x[neigh[3,j]] + x[neigh[4,j]] - 4.0::Float64*x[j]) 
        end
        f[N+1] = ε
        return nothing
end
# Definition of the SDE - stochastic part
function iip_stoc!(f, x, neigh, t)
        for j=1:N
                f[j] = σ
        end
        f[N+1] = 0.0::Float64
        return nothing
end
problem = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T), L.connectivity)

# Compute the forward-time trajectory
solution = solve(problem, EM(), dt=δt)
time = solution.t
states = solution.u

# Dynamic Mode Decomposition
#=
EWS = zeros(Float64, (length(time)-width))
Mode = Array{Float64,2}(undef, N, (length(time)-width))
println("Performing the analysis of the DMD modes")
@showprogress for j=1:(length(time)-width)
        # Assemble the first snapshot matrix
        idx = 1::Int
        X = Array{Float64,2}(undef, N, width)
        for k=j:(j-1)+width
                X[:,idx] = states[k][1:N]
                idx += 1
        end
        # Perform the DMD
        F = svd(X)
        # Assemble the second snapshot matrix
        idx = 1::Int
        Y = Array{Float64,2}(undef, N, width)
        for k=(j+1):j+width
                Y[:,idx] = states[k][1:N]
                idx += 1
        end
        # Perform the eigendecomposition
        S = (F.U)'*Y*(F.Vt)'*inv(diagm(F.S)) 
        Λ = eigvals(S)
        Q = eigvecs(S)
        # Get leading eigenvalue and eigenvector as spatial EWS
        EWS[j] = (findmax(real(Λ)))[1]
        Mode[:,j] = F.U*real(Q[:,(findmax(real(Λ)))[2]])
end
=#

# Covariance analysis
EWS = zeros(Float64, (length(time)-width))
Mode = Array{Float64,2}(undef, N, (length(time)-width))
println("Performing the analysis of the covariance matrix modes")
@showprogress for t=1:(length(time)-width)
        # Get the states of the system across the specified time window 
        state = states[t:((t-1)+width)]
        # Assemble the snasphot matrix for the selected window
        ctr = 1
        X = Array{Float64, 2}(undef, width, N)
        for snapshot in state
                X[ctr, :] = snapshot[1:N]
                ctr +=1
        end
        # Compute the covariance matrix of the snapshots 
        Σ = cov(X, corrected=false)
        #=
        S = Array{Float64, 2}(undef, N, N)
        for j=1:N
                x = X[:,j]
                for k=1:N
                        y = X[:,k]
                        S[j,k] = mean(x.*y)-mean(x)*mean(y) 
                end
        end
        =#
        # Perform the eigendecomposition
        Λ = eigvals(Σ)
        Q = eigvecs(Σ)
        # Get leading eigenvalue and eigenvector as spatial EWS
        EWS[t] = (findmax(real(Λ)))[1]
        Mode[:,t] = real(Q[:,(findmax(real(Λ)))[2]])
end

# Export images of the lattice time evolution
show(L, states, EWS, Mode)
