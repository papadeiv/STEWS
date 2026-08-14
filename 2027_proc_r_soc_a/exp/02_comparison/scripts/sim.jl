"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
Nμ = 100                                                  # Number of parameter values in the sweep
μ_set = collect(LinRange(0.0, 1.20, Nμ))                  # Set of bifurcation parameter values 
ε = 0.0                                                   # Timescale separation
σ = 0.050                                                 # Noise level (additive)
D = (σ^2)/2.0                                             # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ + 4*x - 4*x^3                                # Drift
J(x) = + 4 - 12*x^2                                       # Jacobian 
Λ(t) = ε                                                  # Shift/Ramp
η(x) = σ                                                  # Diffusion

# Simulation parameters
dt = exp10.(range(log10(0.2), log10(0.001), length = 100)) # Timestep size
Nt = round.(Integer, exp10.(range(2, 5, length = 100)))    # Number of timesteps
Ns = length(dt)*length(Nt)                                 # Number of simulations
Ne = 1e2                                                   # Number of particles

# Compute solutions and error approximations of the parameter estimation problems
function generate_samples(stepsize, steps)
        # Initialize stop flag
        stop_flag = false 

        # Parameter sweep loop
        for (idx_μ, μ) in enumerate(μ_set)
                # Compute the relevant equilibria
                equilibria = get_equilibria(f, μ, domain=[-10,10])
                xu = equilibria.unstable[1] 
                xs = equilibria.stable[2]

                # Solve the ensemble problem
                x0 = [xs, μ]
                ensemble = evolve(f, η, Λ, x0, stepsize=stepsize, steps=steps, particles=Ne)

                # Ensemble loop
                data = Matrix{Float64}(undef, convert(Integer, Ne), 3)
                for (idx_sol, solution) in enumerate(ensemble.state)
                        # Truncate and detrend the solution
                        idx_tip = find_tipping(solution, xu)
                        u = solution[1:idx_tip]
                        u = u .- x0[1]

                        # Stop the simulation if a single sample path has tipped
                        if idx_tip < steps 
                                stop_flag = true 
                                @goto tipped 
                        end

                        # Solve the linear least-squares problem of the Euler-Maruyama quasi-likelihood estimation
                        c = (estimate_parameters(u, stepsize))[2]

                        # Compute the lag-1 autocorrelation coefficient ρ and the return rate α
                        α = -100.0
                        ρ = sum(u[1:end-1].*u[2:end])/sum(abs2, u)
                        ρ > 0 && (α = log(ρ)/stepsize)

                        # Compute the sample variance and fit an Ornstein-Uhlenbeck process with linear coefficient θ
                        v = var(u) 
                        θ = -σ^2/(2*v)

                        # Update the data matrix
                        data[idx_sol, :] = [c, α, θ]
                end

                # Export the data
                writeout(data, "$idx_μ/$stepsize/$steps.csv")
        end

        # Return the flag
        @label tipped
        return stop_flag
end
