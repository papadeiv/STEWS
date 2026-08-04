"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")

# Main algorithm 
function main()
        printstyled("Generating and analysing the samples\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for (idx_sim, (stepsize, steps)) in enumerate(Iterators.product(dt, Nt)) 
                # Solve the ensemble problems and export the results
                generate_samples(stepsize, steps)

                # Perform the error analysis for the current simulation setup
                error_analysis(stepsize, steps, idx_sim)
        end

        # Export the error matrix
        writeout(error_matrix, "error_analysis.csv")

        # Plot the error heatmap
        error = reshape(error_matrix[:,3], length(dt), length(Nt))
        logedges(v) = (l = log10.(v); exp10.([l[1] - (l[2]-l[1])/2;
                                      (l[1:end-1] .+ l[2:end]) ./ 2;
                                      l[end] + (l[end]-l[end-1])/2]))

        x, y = logedges(dt), logedges(Nt)
        fig = Figure()
        ax = Axis(fig[1,1], 
                  xscale = log10, yscale = log10,
                  xticks = [1e-3, 1e-2, 1e-1],
                  yticks = [1e+3, 1e+4, 1e+5]
                 )
        hm = heatmap!(ax, x, y, error, colorrange = (minimum(error), maximum(error)), colorscale = log10)
        Colorbar(fig[1,2], hm)

        # Plot the location of the minimizer in the heatmap
        idx_min = argmin(error_matrix[:,3])
        scatter!(ax, error_matrix[idx_min,1], error_matrix[idx_min,2], marker = :star5, color = :yellow, markersize = 20)

        # Export the error heatmap figure
        savefig("error_matrix.png", fig)
end

# Execute the main 
main()
