"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/figs.jl")

# Define the main algorithm
function main()
        # Solve the ensemble problem 
        x0 = [sqrt(-μ0), μ0]
        ensemble = evolve(f, η, Λ, x0, steps=Nt, stepsize=Δt, saveat=dt, particles=Ne)

        # Loop over the ensemble's sample paths
        for (idx_sol, solution) in enumerate(ensemble.state)
                # Plot the timeseries
                N_tip = length(solution)
                lines!(ax1, ensemble.time[1:N_tip], solution[1:N_tip], color = :black, linewidth = 1.0)

                # Identify the tipping and truncate the solution up to the that point
                idx_tip = find_tipping(solution, 0.0)
                t = ensemble.time[1:idx_tip]
                x = solution[1:idx_tip]

                # Compute and plot the quasi-stationary residuals of the truncated timeseries
                residuals = detrend(x; alg = "emd", n_modes = 1).residuals
        end

        # Plot the window-of-oppurtunity
        poly!(ax1, Point2f[(τ,-2), (T,-2), (T,2), (τ,2)], color = (:black,0.10), strokecolor = :black, strokewidth = 2.0)
        poly!(ax2, Point2f[(τ,-1), (T,-1), (T,5), (τ,5)], color = (:black,0.10), strokecolor = :black, strokewidth = 2.0)

        # Format and export the figure
        ax1.limits = (ensemble.time[1], ensemble.time[end], -0.5, +1.2)
        ax2.limits = (ensemble.time[1], ensemble.time[end], -1, +5)
        savefig("extrapolation.pdf", fig)
end

# Execute the main
main()
