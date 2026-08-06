"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/plot.jl")

# Main algorithm 
function main()
        printstyled("Generating and analysing the samples\n"; bold=true, underline=true, color=:light_blue)
        a = [dt[1],dt[end]]
        b = [Nt[1],Nt[end]]
        @showprogress for (idx_sim, (stepsize, steps)) in enumerate(Iterators.product(dt,Nt)) 
                # Solve the ensemble problems and export the results
                generate_samples(stepsize, steps)

                # Perform the statistical analysis for the current simulation setup
                analysis(stepsize, steps, idx_sim)
        end

        # Plot and export the error matrix
        writeout(error_matrix, "error_analysis.csv")
        plot_errormap(error_matrix)
end

# Execute the main 
main()
