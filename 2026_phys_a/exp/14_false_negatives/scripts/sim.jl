"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ0 = -1.0                                                 # Initial bifurcation parameter value
μf = 0.0                                                  # Final bifurcation parameter value
μ1 = -0.5                                                 # Stationary value 1 
μ2 = -0.1                                                 # Stationary value 2
ε = 1e-4                                                  # Timescale separation
σ = 0.050                                                 # Noise level (additive)
D = (σ^2)/2.0                                             # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ - x^2                                        # Drift
Λ(t) = ε                                                  # Shift/Ramp
η(x) = σ                                                  # Diffusion

# Simulation parameters
dt = 1e-1                                                 # Timestep size
Ne = 1e+1                                                 # Number of particles
