"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ_set = collect(LinRange(0.0, 1.0, 100))      # Set of bifurcation parameter values 
ε = 0.0                                       # Timescale separation
σ = 0.200                                     # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ + 4*x - 4*x^3                    # Drift
Λ(t) = ε                                      # Shift/Ramp
η(x) = σ                                      # Diffusion

# Simulation parameters
dt = 5e-2                                     # Timestep
Nt = 1e4                                      # Total number of steps
Ne = 1e1                                      # Number of particles in the ensemble 
