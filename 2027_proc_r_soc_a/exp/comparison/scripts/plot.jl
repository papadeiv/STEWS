# Plot the coefficients' solutions across the parameter sweep 
function plot_solutions(μ_set, outliers_solutions, interquartile_solutions, median_solutions)
        # Create the figure
        fig = Figure()
        ax = Axis(fig[1,1])

        # Plot the error decays
        band!(ax, μ_set, outliers_solutions[:,1], outliers_solutions[:,2], color = (:steelblue, 0.25))
        lines!(ax, μ_set, outliers_solutions[:,1], color = :steelblue, linewidth = 3.0)
        lines!(ax, μ_set, outliers_solutions[:,2], color = :steelblue, linewidth = 3.0)
        lines!(ax, μ_set, interquartile_solutions[:,1], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax, μ_set, interquartile_solutions[:,2], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax, μ_set, median_solutions, color = :red, linewidth = 3.0)

        # Export the figure
        savefig("solutions.png", fig)
end

# Plot spectral decay
function plot_spectral_decay(λ, θm)
        # Create the figure
        fig = Figure(; size = (1200, 400))
        ax1 = Axis(fig[1,1])
        ax2 = Axis(fig[1,2])
        ax3 = Axis(fig[1,3])

        # Plot the spectrum 
        scatter!(ax1, μ_set, λ[:,1], markersize = 15, color = (:red,0.5)) 
        scatter!(ax2, μ_set, λ[:,2], markersize = 15, color = (:red,0.5)) 
        scatter!(ax3, μ_set, λ[:,3], markersize = 15, color = (:red,0.5))

        # Export the figure
        savefig("spectrum.png", fig)
end

# Plot the error decays with the bifurcation parameter 
function plot_error_decay(μ_set, outliers_error, interquartile_error, median_error)
        # Create the figure
        fig = Figure()
        ax = Axis(fig[1,1], yscale = log10)

        # Plot the error decays
        band!(ax, μ_set, outliers_error[:,1], outliers_error[:,2], color = (:steelblue, 0.25))
        lines!(ax, μ_set, outliers_error[:,1], color = :steelblue, linewidth = 3.0)
        lines!(ax, μ_set, outliers_error[:,2], color = :steelblue, linewidth = 3.0)
        lines!(ax, μ_set, interquartile_error[:,1], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax, μ_set, interquartile_error[:,2], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax, μ_set, median_error, color = :red, linewidth = 3.0)

        # Export the figure
        savefig("error.png", fig)
end
