"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ0 = -0.5                                         # Initial value of the bifurcation parameter
ε = 1e-4                                          # Timescale separation
σ = 0.050                                         # Noise level (additive)
D = (σ^2)/2.0                                     # Diffusion level (additive) 
τ = -μ0/ε                                         # Breaking time
T = -μ0/ε + 1.018793*ε^(-1/3)                     # Tipping time (point of no-return)
T1, T2 = 0.25, 0.75                               # Points for the inference window

# Dynamical system  
f(x, μ) = -μ - x^2                                # Drift
J(x) = -2*x                                       # Jacobian 
Λ(t) = ε                                          # Shift/Ramp
η(x) = σ                                          # Diffusion

# Simulation parameters
dt = 1e-1                                         # Timestep
Nt = 5e+4                                         # Number of steps 
Ne = 1e+2                                         # Number of particles in the ensemble 

# Solve the dynamic saddle problems and export the results 
function generate_samples()
        # Solve the ensemble problem 
        x0 = [sqrt(-μ0), μ0]
        ensemble = evolve(f, η, Λ, x0, steps=Nt, stepsize=dt, particles=Ne)

        # Determine the earliest tipping in the ensemble 
        earliest_tip_idx = length(ensemble.time) 
        for (idx_sol, solution) in enumerate(ensemble.state)
                # Plot the timeseries 
                idx_sol == 1 && lines!(ax1, ensemble.time[1:length(solution)], solution, color = :black, linewidth = 1.0)

                # Extract the tipping time
                idx_tip = find_tipping(solution, 0.0)
                idx_tip < earliest_tip_idx && (earliest_tip_idx = idx_tip) 
        end

        # Compute the inference window time index
        ΔT = [convert(Integer, floor(T1*earliest_tip_idx)), convert(Integer, floor(T2*earliest_tip_idx))]

        # Debug the simulation
        println("Longest simulation time:       $(ensemble.time[end]) (idx: $(length(ensemble.time)))")
        println("Longest simulation parameter:  $(ensemble.parameter[end]) (idx: $(length(ensemble.time)))")
        println("Earliest tipping time:         $(ensemble.time[earliest_tip_idx]) (idx: $(earliest_tip_idx))")
        println("Earliest tipping parameter:    $(ensemble.parameter[earliest_tip_idx]) (idx: $(earliest_tip_idx))")
        println("Inference window time:         [$(ensemble.time[ΔT[1]]), $(ensemble.time[ΔT[2]])] (length: $(length(ensemble.time[ΔT[1]:ΔT[2]])))")
        println("Inference window parameter:    [$(ensemble.parameter[ΔT[1]]), $(ensemble.parameter[ΔT[2]])] (length: $(length(ensemble.parameter[ΔT[1]:ΔT[2]])))")
        println("Sliding window length:         $(ΔT[1])")

        return (
                ensemble = ensemble,
                tipping = earliest_tip_idx,
                window = ΔT
               )
end
