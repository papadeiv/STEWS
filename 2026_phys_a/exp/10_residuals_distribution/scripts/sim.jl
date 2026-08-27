"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ0 = 0.0                                                 # Initial value of the bifurcation parameter
μf = 0.2                                                 # Final value of the bifurcation parameter
ε = 5e-4                                                  # Timescale separation
σ = 0.100                                                 # Noise level (additive)
D = (σ^2)/2.0                                             # Diffusion level (additive) 
dt = 1e-1                                                 # Timestep size
Ne = 1e+3                                                 # Number of particles

# Dynamical system  
f(x, μ) = -μ + 4*x - 4*x^3                                # Drift
Λ(t) = ε                                                  # Shift/Ramp
η(x) = σ                                                  # Diffusion
