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
        # Solve the nonautonomous problems
        samples = generate_samples()
        ensemble = samples.ensemble
        ΔT = samples.window

        # Estimate the return rate in the inference window 
        Ns = 0 
        for (idx_sol, solution) in enumerate(ensemble.state)
                # Extract and truncate the timestamps and the timeseries in the inference window
                t = ensemble.time[1:ΔT[2]]
                x = solution[1:ΔT[2]]

                # Compute the quasi-stationary residuals of the truncated timeseries
                x = detrend(x; alg = "linear", timestamps = t).residuals

                # Estimate the return rate EWS across a sliding window
                Ns = analyse(t, x, ΔT[1], idx_sol)
        end

        # Extrapolate the tipping time (zero-crossing of the return rate in the future)
        extrapolate(ensemble.time, Ns)

        # Plot the inference window and the tipping threshold
        poly!(ax1, Point2f[(ensemble.time[ΔT[1]],-0.5), (ensemble.time[ΔT[2]],-0.5), (ensemble.time[ΔT[2]],1.2), (ensemble.time[ΔT[1]],1.2)], color = (:black,0.10), strokecolor = :black, strokewidth = 2.0)
        poly!(ax2, Point2f[(ensemble.time[ΔT[1]],-0.5), (ensemble.time[ΔT[2]],-0.5), (ensemble.time[ΔT[2]],7), (ensemble.time[ΔT[1]],7)], color = (:black,0.10), strokecolor = :black, strokewidth = 2.0)
        poly!(ax3, Point2f[(ensemble.time[ΔT[1]],-0.5), (ensemble.time[ΔT[2]],-0.5), (ensemble.time[ΔT[2]],7), (ensemble.time[ΔT[1]],7)], color = (:black,0.10), strokecolor = :black, strokewidth = 2.0)
        poly!(ax4, Point2f[(ensemble.time[ΔT[1]],-0.5), (ensemble.time[ΔT[2]],-0.5), (ensemble.time[ΔT[2]],7), (ensemble.time[ΔT[1]],7)], color = (:black,0.10), strokecolor = :black, strokewidth = 2.0)
        vlines!(ax1, T, color = :black, linestyle = :dash, linewidth = 2.0)
        vlines!(ax2, T, color = :black, linestyle = :dash, linewidth = 2.0)
        vlines!(ax3, T, color = :black, linestyle = :dash, linewidth = 2.0)
        vlines!(ax4, T, color = :black, linestyle = :dash, linewidth = 2.0)

        # Format and export the figure
        ax1.limits = (ensemble.time[1], ensemble.time[end], -0.5, +1.2)
        ax2.limits = (ensemble.time[1], ensemble.time[end], -0.5, 7)
        ax3.limits = (ensemble.time[1], ensemble.time[end], -0.5, 7)
        ax4.limits = (ensemble.time[1], ensemble.time[end], -0.5, 7)
        savefig("extrapolation.pdf", fig)
end

# Execute the main
main()
