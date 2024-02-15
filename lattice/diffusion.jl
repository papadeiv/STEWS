using DifferentialEquations, LaTeXStrings, ProgressMeter
include("lattice.jl")

# Define the lattice
n = 15
m = 20
N = n*m
L = Lattice(n, m)

# Initial condition of the lattice
x0 = 0.0::Float64.*ones(Float64, length(L.grid))
x0[50] = 1.0::Float64
# Initial condition of the bifurcation parameter (degree of smooth muscle activation)
push!(x0, 0.0::Float64)

# Model's parameters
σ = 0.50::Float64 # noise level
ε = 1.00::Float64 # parameter's drift
k = 14.1::Float64 # (airway) smooth muscle mass
A = 0.63::Float64 # inter-airway coupling
Pi = 0.96::Float64 # inflection point of pressure-radius interdependence
Pb = 0.0::Float64 # breathing pressure
Pb0 = 7.25::Float64 # breating pressure's IC

# Sliding window's parameters
T = 2.000
δt = 1e-3

function iip_det!(f, x, neigh, t)
        sum = 0.0::Float64
        for j=1:N
                sum += x[j]^4
        end
        Pb = (Pb0*N)/(sum)
        println("\n")
        #println(sum)
        for j=1:N
                r = x[j]
                No = x[neigh[1,j]]
                W = x[neigh[2,j]]
                S = x[neigh[3,j]]
                E = x[neigh[4,j]] 
                #println("$(j):  x = $(r), N($(neigh[1,j])) = $(No), E($(neigh[4,j])) = $(E), S($(neigh[3,j])) = $(S), W($(neigh[2,j])) = $(W)")
                f[j] = x[N+1]*(No + W + S + E)
        end
        f[N+1] = ε
        return nothing
end
function iip_stoc!(f, x, neigh, t)
        for j=1:N
                f[j] = σ
        end
        f[N+1] = 0.0::Float64
        return nothing
end
problem = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T), L.connectivity)
solution = solve(problem, EM(), dt=δt)

time = solution.t
states = solution.u

show(L, states)
