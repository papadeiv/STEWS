using DifferentialEquations, LaTeXStrings, ProgressMeter
include("Lattice.jl")

# Define the lattice
Nx = 10 
Ny = 10 
N = Nx*Ny
L = Lattice(Nx, Ny)

# Initial condition of the lattice
x0 = ones(Float64, length(L.grid))
# Initial condition of the bifurcation parameter (degree of smooth muscle activation)
push!(x0, 1.0::Float64)

# Model's parameters
σ = 0.1::Float64
k = 14.1::Float64 # (airway) smooth muscle mass
A = 0.63::Float64 # inter-airway coupling
Pi = 0.96::Float64 # inflection point of pressure-radius interdependence
Pb = 0.0::Float64 # breathing pressure
Pb0 = 7.25::Float64 # breating pressure's IC

# Sliding window's parameters
T = 3.00
δt = 1e-3

function iip_det!(f, x, neigh, t)
        sum = 0.0::Float64
        for j=1:N
                sum += x[j]^4
        end
        Pb = (Pb0*N)/(sum)
        for j=1:N
                f[j] = (1+ exp(-Pb + x[end]*(k/x[j]) + Pb*A*(x[j]^4 + x[neigh[1,j]]^4 + x[neigh[2,j]]^4 + x[neigh[3,j]]^4 + x[neigh[4,j]]^4)*(1.0::Float64 - x[j] + 1.5::Float64*(1-x[j])^2) + Pi))^(-1) - x[j]
        end
        f[end] = 1.0::Float64
        return nothing
end
function iip_stoc!(f, x, neigh, t)
        for j=1:N
                f[j] =  σ
        end
        f[end] = 0.0::Float64
        return nothing
end
problem = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T), L.connectivity)
solution = solve(problem, EM(), dt=δt)

time = solution.t
states = solution.u

show(L, states)
