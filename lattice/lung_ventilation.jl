using DifferentialEquations, LaTeXStrings, ProgressMeter
include("lattice.jl")

# Define the lattice
n = 3
m = 3
N = n*m
L = Lattice(n, m)

# Initial condition of the lattice
x0 = 1.0::Float64.*ones(Float64, length(L.grid))
# Initial condition of the bifurcation parameter (degree of smooth muscle activation)
push!(x0, 0.0::Float64)
# SDE's parameters
σ = 0.05::Float64 # noise level
push!(x0, σ)
ε = 0.10::Float64 # parameter's drift
push!(x0, ε)

# Model's parameters
k = 14.1::Float64 # (airway) smooth muscle mass
A = 0.63::Float64 # inter-airway coupling
Pi = 0.96::Float64 # inflection point of pressure-radius interdependence
Pb = 0.0::Float64 # breathing pressure
Pb0 = 7.25::Float64 # breating pressure's IC

# Sliding window's parameters
T = 01.00
δt = 1e-2
width = 20

# Definition of the dynamics (lung ventilation) function
sigmoid(x::Real) = one(x)/(one(x) + exp(-x))
# Definition of the SDE - deterministic part
function iip_det!(f, x, neigh, t)
        sum = 0.0::Float64
        for j=1:N
                sum += x[j]^4
        end
        Pb = (Pb0*N)/(sum)
        for j=1:N
                r = -Pb + x[N+1]*(k/x[j]) - Pb*A*(x[j]^4 + x[neigh[1,j]]^4 + x[neigh[2,j]]^4 + x[neigh[3,j]]^4 + x[neigh[4,j]]^4)*(1.0::Float64 - x[j] + 1.5::Float64*(1-x[j])^2) + Pi
                f[j] = sigmoid(-r) - x[j]
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
for j=1:(length(time)-width)
        # Assemble snapshot matrices
        idx = 1::Int
        X = Array{Float64,2}(undef, N, width)
        for k=j:(j-1)+width
                X[:,idx] = states[k][1:N]
                idx += 1
        end
        println(j)
        println(X)
        println("\n")
end

# Export images of the lattice time evolution
#show(L, states)
