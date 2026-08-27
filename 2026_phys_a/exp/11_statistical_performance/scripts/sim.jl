"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# Simulation parameters
Nh = 150                                                  # Number of grid points in the simulation
Ns = Nh^2                                                 # Number of simulations
dt = exp10.(range(log10(0.1), log10(0.01), length = Nh))  # Timestep size
Nt = round.(Integer, exp10.(range(3, 4, length = Nh)))    # Number of timesteps
Ne = 1e3                                                  # Number of particles

# System parameters
μ_set = [1.1, 1.2, 1.3]                                   # Set of bifurcation parameter values
ε = exp10.(range(log10(1e-5), log10(1e-3), length = Nh))  # Timescale separation
σ = exp10.(range(log10(1e-3), log10(1e-1), length = Nh))  # Noise level (additive)
f(x, μ) = -μ + 4*x - 4*x^3                                # Drift

# Simulation's parameters loop
function generate_samples_sim(stepsize, steps)
        # Initialize the dynamics
        Λ(t) = 0.000
        η(x) = 0.050

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
                data = Vector{Float64}(undef, convert(Integer, Ne))
                for (idx_sol, solution) in enumerate(ensemble.state)
                        # Detrend the solution
                        u = solution .- x0[1]

                        # Solve the LLS problem and compute the modified escape rate 
                        θ = solve_lls(u, stepsize)
                        data[idx_sol] = compute_ews(θ)
                end

                # Export the data
                writeout(data, "sim_par/$idx_μ/$stepsize/$steps.csv")
        end
end

# System's parameters loop
function generate_samples_sys(noise, speed)
        # Initialize the dynamics
        Λ(t) = speed 
        η(x) = noise

        # Initialize the simulation parameters
        local Nt = 1e+4
        local dt = 1e-1

        # Parameter sweep loop
        for (idx_μ, μf) in enumerate(μ_set)
                # Compute the starting parameter and the initial condition
                μ0 = μf - speed*(Nt*dt) 
                x0 = [(get_equilibria(f, μ0, domain=[-10,10])).stable[2], μ0]


                # Solve the ensemble problem
                ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)
                t = ensemble.time
                μ = ensemble.parameter

                # Ensemble loop
                data = Vector{Float64}(undef, convert(Integer, Ne))
                for (idx_sol, solution) in enumerate(ensemble.state)
                        # Detrend the solution
                        u = detrend(solution; alg = "linear", timestamps = t).residuals

                        # Solve the LLS problem and compute the modified escape rate 
                        θ = solve_lls(u, 1e-1)
                        data[idx_sol] = compute_ews(θ)
                end

                # Export the data
                writeout(data, "sys_par/$idx_μ/$noise/$speed.csv")
        end
end
