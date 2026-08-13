# Plot the coefficients' solutions across the parameter sweep 
function plot_solutions(μ_set, outliers_solutions, interquartile_solutions, median_solutions, outliers_eigenvalue, interquartile_eigenvalue, median_eigenvalue)
        # Create the figure
        fig = Figure(; size = (600,300))
        ax1 = Axis(fig[1,1], xlabel = L"\mu", ylabel = "inv. leading eigenvalue", title = "From LLS solutions")#, limits = (μ_set[1],μ_set[end],0,6))
        ax2 = Axis(fig[1,2], xlabel = L"\mu", title = "From AR(1) fit")#, limits = (new_μ_set[1],new_μ_set[end],0,6))

        # Plot the leading eigenvalue reconstruction from LLS 
        band!(ax1, μ_set, outliers_solutions[2][:,1], outliers_solutions[2][:,2], color = (:steelblue, 0.25))
        lines!(ax1, μ_set, outliers_solutions[2][:,1], color = :steelblue, linewidth = 3.0)
        lines!(ax1, μ_set, outliers_solutions[2][:,2], color = :steelblue, linewidth = 3.0)
        lines!(ax1, μ_set, interquartile_solutions[2][:,1], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax1, μ_set, interquartile_solutions[2][:,2], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax1, μ_set, median_solutions[2], color = :red, linewidth = 3.0)

        # Plot the leading eigenvalue reconstruction from AR(1) 
        band!(ax2, μ_set, outliers_eigenvalue[:,1], outliers_eigenvalue[:,2], color = (:steelblue, 0.25))
        lines!(ax2, μ_set, outliers_eigenvalue[:,1], color = :steelblue, linewidth = 3.0)
        lines!(ax2, μ_set, outliers_eigenvalue[:,2], color = :steelblue, linewidth = 3.0)
        lines!(ax2, μ_set, interquartile_eigenvalue[:,1], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax2, μ_set, interquartile_eigenvalue[:,2], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax2, μ_set, median_eigenvalue, color = :red, linewidth = 3.0)

        # Loop over the parameter values
        new_μ_set = LinRange(0.0, 1.5396, Nμ)
        λ = Vector{Real}(undef, Nμ)
        for (idx_μ, μ) in enumerate(new_μ_set)
                # Compute the stable equilibrium at the current parameter value
                x0 = (get_equilibria(f, μ, domain=[-10,10])).stable[2]

                # Evaluate the Jacobian at the current equilibrium
                λ[idx_μ] = J(x0)
        end
        
        # Plot the ground truth of the leading eigenvalue 
        lines!(ax1, new_μ_set, λ, color = :black, linewidth = 2.0)
        lines!(ax2, new_μ_set, λ, color = :black, linewidth = 2.0)

        # Export the figure
        savefig("comparison.png", fig)
end
