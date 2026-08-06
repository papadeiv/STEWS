"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Ground truth of the Langevin dynamics and its Taylor approximation
U(x, μ) = 1 + μ*x - 2*(x^2) + x^4
Ux(x, μ) = μ -4*x + 4*x^3
U2x(x, μ) = -4 + 12*x^2
U3x(x, μ) = + 24*x
T(x, x0, μ) = U(x0, μ) + Ux(x0, μ)*(x - x0) + (U2x(x0, μ)/2)*(x - x0)^2 + (U3x(x0, μ)/6)*(x - x0)^3

# Error data arrays
median_error = Vector{Real}(undef, convert(Integer, Nμ))
interquartile_error = Matrix{Real}(undef, convert(Integer, Nμ), 2)
outliers_error = Matrix{Real}(undef, convert(Integer, Nμ), 2)
error_matrix = Matrix{Real}(undef, length(dt)*length(Nt), 3)

# Perform the statistical analysis of the generated samples
function analysis(stepsize, steps, idx_sim)
        # Initialize data structures for the analysis
        ellipsoids = Vector{NamedTuple}(undef, Nμ) 
        eigenbases = Vector{Matrix{Real}}(undef, Nμ)
        centroids = Matrix{Real}(undef, Nμ, 3)
        natural_spectrum = Matrix{Real}(undef, Nμ, 3)
        centered_spectrum = Matrix{Real}(undef, Nμ, 3)
        mean_sol_norm = Vector{Real}(undef, Nμ)

        # Loop over the parameter values
        for (idx_μ, μ) in enumerate(μ_set)
                # Import the data
                data = readin("$stepsize/$steps/$idx_μ.csv")
                θ = data[:,1:3]

                # Compute the data ellipsoid
                ellipsoid = fit_ellipsoid(θ)
                ellipsoids[idx_μ] = ellipsoid.ellipsoid
                eigenbases[idx_μ] = ellipsoid.eigenbasis
                centroids[idx_μ,:] = ellipsoid.mean
                mean_sol_norm[idx_μ] = sum((ellipsoid.mean).^2)

                # Compute the spectrum of the covariance of the sample
                centered_spectrum[idx_μ,:] = ellipsoid.spectrum 
                #θ = θ .- mean(θ; dims=1)
                Σ = θ'*θ/(size(θ,1) - 1)
                Λ = eigen(Symmetric(Matrix(Σ)))
                λ = Λ.values
                v = Λ.vectors
                natural_spectrum[idx_μ,:] = λ 

                # Compute the error statistics
                error = data[:,4]
                median_error[idx_μ] = median(error)
                interquartile_error[idx_μ,:] = quantile(error, [0.25, 0.75])
                outliers_error[idx_μ,:] = collect(extrema(error))
        end

        # Update the error matrix
        error_matrix[idx_sim,:] = [stepsize; steps; median_error[idx_analysis]]

        # Plot the data ellipsoids
        if stepsize == 0.1 && steps == 100000
                plot_ellipsoids(ellipsoids, eigenbases, centroids, μ_set, idx_sim, stepsize, steps)
                plot_error_decay(μ_set, outliers_error, interquartile_error, median_error, idx_sim, stepsize, steps)
        end
        
        # Plot and export the solutions, error and spectral decay figures
        plot_mean_solutions(centroids, μ_set, idx_sim, stepsize, steps)
        plot_spectral_decay(natural_spectrum, centered_spectrum, mean_sol_norm, μ_set, idx_sim, stepsize, steps)
        #plot_error_decay(μ_set, outliers_error, interquartile_error, median_error, idx_sim, stepsize, steps)
end
