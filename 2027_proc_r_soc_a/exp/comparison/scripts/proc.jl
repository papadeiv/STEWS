"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Flag for tipped states
tipped = false

# Ground truth of the Langevin dynamics and its Taylor approximation
U(x, μ) = 1 + μ*x - 2*(x^2) + x^4
Ux(x, μ) = μ -4*x + 4*x^3
U2x(x, μ) = -4 + 12*x^2
U3x(x, μ) = + 24*x
T(x, x0, μ) = U(x0, μ) + Ux(x0, μ)*(x - x0) + (U2x(x0, μ)/2)*(x - x0)^2 + (U3x(x0, μ)/6)*(x - x0)^3

# Initialize data structures for the analysis
median_error = Vector{Real}(undef, convert(Integer, Nμ))
interquartile_error = Matrix{Real}(undef, convert(Integer, Nμ), 2)
outliers_error = Matrix{Real}(undef, convert(Integer, Nμ), 2)
median_solutions = Vector{Real}(undef, convert(Integer, Nμ))
interquartile_solutions = Matrix{Real}(undef, convert(Integer, Nμ), 2)
outliers_solutions = Matrix{Real}(undef, convert(Integer, Nμ), 2)
spectrum = Matrix{Real}(undef, Nμ, 3)

# Perform the statistical analysis of the generated samples
function analysis()
        # Loop over the parameter values
        printstyled("Analyzing the samples\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for (idx_μ, μ) in enumerate(μ_set)
                # Import the data
                data = readin("solutions/$idx_μ.csv")
                θ = data[:,1:3]

                # Compute the solutions statistics
                θ2 = θ[:,2]
                median_solutions[idx_μ] = median(θ2)
                interquartile_solutions[idx_μ,:] = quantile(θ2, [0.25, 0.75])
                outliers_solutions[idx_μ,:] = collect(extrema(θ2))

                # Compute the data ellipsoid
                ellipsoid = fit_ellipsoid(θ)
                spectrum[idx_μ,:] = ellipsoid.spectrum 

                # Compute the error statistics
                error = data[:,4]
                median_error[idx_μ] = median(error)
                interquartile_error[idx_μ,:] = quantile(error, [0.25, 0.75])
                outliers_error[idx_μ,:] = collect(extrema(error))
        end

        # Plot and export the solutions, error and spectral decay figures
        plot_solutions(μ_set, outliers_solutions, interquartile_solutions, median_solutions)
        plot_spectral_decay(spectrum, μ_set)
        plot_error_decay(μ_set, outliers_error, interquartile_error, median_error)
end
