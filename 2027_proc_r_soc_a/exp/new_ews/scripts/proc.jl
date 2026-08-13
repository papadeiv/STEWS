"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Flag for tipped states
tipped = false

# Initialize data structures for the analysis
median_eigenvalue = Vector{Real}(undef, Nμ)
interquartile_eigenvalue = Matrix{Real}(undef, Nμ, 2)
outliers_eigenvalue = Matrix{Real}(undef, Nμ, 2)
median_solutions = [zeros(Nμ) for _ in 1:3]
interquartile_solutions = [zeros(Nμ, 2) for _ in 1:3]
outliers_solutions = [zeros(Nμ, 2) for _ in 1:3]

# Perform the statistical analysis of the generated samples
function analysis()
        # Loop over the parameter values
        printstyled("Analyzing the samples\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for (idx_μ, μ) in enumerate(μ_set)
                # Import the data 
                data = readin("solutions/$idx_μ.csv")
                θ = data[:,1:3]
                λ = data[:,4]

                # Perform the coefficients analysis
                for (idx_col, col) in enumerate(eachcol(θ))
                        median_solutions[idx_col][idx_μ] = median(col)
                        interquartile_solutions[idx_col][idx_μ,:] = quantile(col, [0.25, 0.75])
                        outliers_solutions[idx_col][idx_μ,:] = collect(extrema(col))
                end

                # Perform the AR(1) regression analysis
                median_eigenvalue[idx_μ] = median(λ)
                interquartile_eigenvalue[idx_μ,:] = quantile(λ, [0.25, 0.75])
                outliers_eigenvalue[idx_μ,:] = collect(extrema(λ))
        end

        # Plot and export the solutions, error and spectral decay figures
        plot_solutions(μ_set, outliers_solutions, interquartile_solutions, median_solutions, outliers_eigenvalue, interquartile_eigenvalue, median_eigenvalue)
end
