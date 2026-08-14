"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Flag for tipped states
tipped = false

# Initialize data structures for the analysis
steps_minimizer = Matrix{Real}(undef, Nμ, 4)
stepsize_minimizer = Matrix{Real}(undef, Nμ, 4)

# Perform the statistical analysis of the generated samples
function analysis()
        # Parameter sweep loop
        printstyled("Analyzing the samples\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for (idx_μ, μ) in enumerate(μ_set)
                # Compute the ground truth of the leading eigenvalue
                xs = (get_equilibria(f, μ, domain=[-10,10])).stable[2]
                λ = J(xs)

                # Simulation parameters loop
                relative_errors = Matrix{Real}(undef, Ns, 3)
                for (idx_sim, (stepsize, steps)) in enumerate(Iterators.product(dt,Nt))
                        # Import the data
                        data = readin("$idx_μ/$stepsize/$steps.csv")

                        # Compute the median of the three estimates
                        c = median(data[:,1])
                        α = median(data[:,2])
                        θ = median(data[:,3])

                        # Compute the relative errors between the median estimates and the ground truth
                        relative_errors[idx_sim,1] = abs(λ-c)/abs(λ)
                        relative_errors[idx_sim,2] = abs(λ-α)/abs(λ)
                        relative_errors[idx_sim,3] = abs(λ-θ)/abs(λ)
                end

                # Extract the minimum and maximum relative errors
                err_min, err_max = minimum(relative_errors), maximum(relative_errors)

                # Create the figure
                include("./scripts/figs.jl")

                # Transformation of the coordinate axes in logscale so that it does not interpolate negative values
                logedges(v) = (l = log10.(v); exp10.([l[1] - (l[2]-l[1])/2;
                              (l[1:end-1] .+ l[2:end]) ./ 2;
                              l[end] + (l[end]-l[end-1])/2]))
                x, y = logedges(dt), logedges(Nt)

                # Plot the heatmap of the relative errors and the minimizers
                stepsize_minimizer[idx_μ,1] = μ 
                steps_minimizer[idx_μ,1] = μ 
                
                # LLS -----------------------------------------------------------------------------------------------------
                error_map = reshape(relative_errors[:,1], length(dt), length(Nt))
                heatmap!(ax1, x, y, error_map, colorrange = (err_min, err_max))#, colorscale = log10)

                idx_min = argmin(error_map)
                scatter!(ax1, dt[idx_min[1]], Nt[idx_min[2]], marker = :star5, color = :yellow, markersize = 20)

                stepsize_minimizer[idx_μ,2] = dt[idx_min[1]]
                steps_minimizer[idx_μ,2] = Nt[idx_min[2]]

                # AR(1) ---------------------------------------------------------------------------------------------------
                error_map = reshape(relative_errors[:,2], length(dt), length(Nt))
                heatmap!(ax2, logedges(dt), logedges(Nt), error_map, colorrange = (err_min, err_max))#, colorscale = log10)

                idx_min = argmin(error_map)
                scatter!(ax2, dt[idx_min[1]], Nt[idx_min[2]], marker = :star5, color = :yellow, markersize = 20)

                stepsize_minimizer[idx_μ,3] = dt[idx_min[1]]
                steps_minimizer[idx_μ,3] = Nt[idx_min[2]]

                # OUP -----------------------------------------------------------------------------------------------------
                error_map = reshape(relative_errors[:,3], length(dt), length(Nt))
                heatmap!(ax3, logedges(dt), logedges(Nt), error_map, colorrange = (err_min, err_max))#, colorscale = log10)

                idx_min = argmin(error_map)
                scatter!(ax3, dt[idx_min[1]], Nt[idx_min[2]], marker = :star5, color = :yellow, markersize = 20)

                stepsize_minimizer[idx_μ,4] = dt[idx_min[1]]
                steps_minimizer[idx_μ,4] = Nt[idx_min[2]]

                # Plot the location of the setup of Figure 1 in the error map
                scatter!(ax1, 1e-1, 1e+3, color = :brown1, markersize = 20)
                scatter!(ax2, 1e-1, 1e+3, color = :brown1, markersize = 20)
                scatter!(ax3, 1e-1, 1e+3, color = :brown1, markersize = 20)

                # Add the colorbar 
                Colorbar(fig[1,4], 
                         size = 25,
                         label = "relative error",
                         limits = (trunc(err_min,digits=1), trunc(err_max,digits=1)), 
                         ticks = [trunc(err_min,digits=1), trunc(err_max,digits=1)], 
                        )
                # Export the figure
                savefig("02_comparison/$idx_μ.pdf", fig)
        end

        # Export the minimizers
        writeout(stepsize_minimizer, "stepsize_minimizer.csv")
        writeout(steps_minimizer, "steps_minimizer.csv")
end
