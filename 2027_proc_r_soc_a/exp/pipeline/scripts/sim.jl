"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
Nμ = 20                                                   # Number of parameter values in the sweep
μ_set = collect(LinRange(0.0, 1.0, Nμ))                   # Set of bifurcation parameter values 
idx_analysis = findfirst(≥(0.5), μ_set)                   # Parameter value for the error analysis 
ε = 0.0                                                   # Timescale separation
σ = 0.200                                                 # Noise level (additive)
D = (σ^2)/2.0                                             # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ + 4*x - 4*x^3                                # Drift
Λ(t) = ε                                                  # Shift/Ramp
η(x) = σ                                                  # Diffusion

# Simulation parameters
Ne = 1e2                                                   # Number of particles in the ensemble 
dt = exp10.(range(log10(0.1), log10(0.001), length = 20))  # Timestep size
#dt = [1e-1,1e-2]
Nt = round.(Integer, exp10.(range(3, 5, length = 20)))     # Number of timesteps
#Nt = [1e+4,1e+5]

# Compute solutions and error approximations of the parameter estimation problems
function generate_samples(stepsize, steps)
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
                data = Matrix{Float64}(undef, convert(Integer, Ne), 4)
                for (idx_sol, solution) in enumerate(ensemble.state)
                        # Truncate and detrend the solution
                        idx_tip = find_tipping(solution, xu)
                        u = solution[1:idx_tip]
                        u = u .- x0[1]

                        # Linear least-squares solution of the Euler-Maruyama quasi-likelihood estimation 
                        θ = estimate_parameters(u, stepsize)
                        try 
                                V = fit_potential(θ, xs, U, μ)
                                error = get_error(U, V, (xs, xu), μ)
                                data[idx_sol, :] = [θ; error]
                        catch
                                data[idx_sol, :] = [θ; Inf]
                        end
                end

                # Clean and export the data
                E∞ = maximum(filter(isfinite, data[:,4]))
                for idx_err in eachindex(data[:,4])
                        isinf(data[idx_err,4]) && (data[idx_err,4] = E∞)
                end
                writeout(data, "$stepsize/$steps/$idx_μ.csv")
        end
end
