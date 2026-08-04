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

# Perform error analysis on the generated samples
function error_analysis(stepsize, steps, idx_sim)

        fig = Figure()
        ax = Axis3(fig[1,1], azimuth = 0.2π, aspect = (2, 3, 1), 
                   xlabel = L"\theta_1",
                   ylabel = L"\theta_2",
                   zlabel = L"\theta_3",
                   title = "dt=$stepsize, Nt=$steps")

        # Loop over the parameter values
        for (idx_μ, μ) in enumerate(μ_set)
                # Import the data
                data = readin("$stepsize/$steps/$idx_μ.csv")
                θ = data[:,1:3]
                ellipsoid = fit_ellipsoid(θ)

                # Plot the uncertainty ellipsoid 
                if idx_μ == 1 || idx_μ == Nμ 
                        surface!(ax, ellipsoid.ellipsoid.X, ellipsoid.ellipsoid.Y, ellipsoid.ellipsoid.Z, 
                         color = idx_μ, colormap = :berlin, colorrange = (1, Nμ))
                else
                        surface!(ax, ellipsoid.ellipsoid.X, ellipsoid.ellipsoid.Y, ellipsoid.ellipsoid.Z, 
                         color = idx_μ, colormap = (:berlin, 0.01), colorrange = (1, Nμ), rasterize = 4)
                end

                # Update the error statistics
                error = data[:,4]
                median_error[idx_μ] = median(error)
                interquartile_error[idx_μ,:] = quantile(error, [0.25, 0.75])
                outliers_error[idx_μ,:] = collect(extrema(error))
        end

        # Export the ellipsoid 3D figure
        Colorbar(fig[2,1], limits = (μ_set[1], μ_set[end]), colormap = :berlin, vertical = false, label = L"\mu")
        savefig("ellipsoid/$idx_sim.png", fig)

        # Update the error matrix
        error_matrix[idx_sim,:] = [stepsize; steps; median_error[idx_analysis]]

        # Plot the error vs parameter value graph
        fig = Figure()
        ax = Axis(fig[1,1], yscale = log10, title = "dt=$stepsize, Nt=$steps")
        band!(ax, μ_set, outliers_error[:,1], outliers_error[:,2], color = (:steelblue, 0.25))
        lines!(ax, μ_set, outliers_error[:,1], color = :steelblue, linewidth = 3.0)
        lines!(ax, μ_set, outliers_error[:,2], color = :steelblue, linewidth = 3.0)
        lines!(ax, μ_set, interquartile_error[:,1], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax, μ_set, interquartile_error[:,2], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax, μ_set, median_error, color = :red, linewidth = 3.0)
        savefig("error_decay/$idx_sim.png", fig)
end
