# Plot the data ellipsoid 
function plot_ellipsoids(ellipsoids, eigenbases, centroids, μ_set, idx_sim, stepsize, steps)
        # Create the figure
        fig = Figure()#; size = (1200, 1200))
        ax = Axis3(fig[1,1], azimuth = 0.2π, aspect = (2, 3, 1), 
                   title = "dt=$stepsize, Nt=$steps",
                   xlabel = L"\theta_1",
                   ylabel = L"\theta_2",
                   zlabel = L"\theta_3",
                  )
        Colorbar(fig[2,1], limits = (μ_set[1], μ_set[end]), colormap = :berlin, vertical = false, label = L"\mu")

        # Loop over the ellipsoids
        for (idx_ell, ellipsoid) in enumerate(ellipsoids)
                # Check which step in the parameter sweep you are
                if idx_ell == 1 
                        # Extract the center of mass of the ellipsoid
                        centroid = centroids[idx_ell,:]
                        # Plot the eigenbasis 
                        for eigenvector in eachcol(eigenbases[idx_ell])
                                lines!(ax, [Point3f(centroid), Point3f(centroid .+ 4.0.*eigenvector)]; color = :red, linewidth = 3)
                        end
                        # Plot the starting ellipsoid fully opaque
                        surface!(ax, ellipsoid.X, ellipsoid.Y, ellipsoid.Z, 
                                 color = idx_ell, colormap = (:berlin,0.5), colorrange = (1, Nμ))
                elseif idx_ell == Nμ
                        # Extract the center of mass of the ellipsoid
                        centroid = centroids[idx_ell,:]
                        # Plot the eigenbasis 
                        for eigenvector in eachcol(eigenbases[idx_ell])
                                lines!(ax, [Point3f(centroid), Point3f(centroid .+ 4.0.*eigenvector)]; color = :red, linewidth = 3)
                        end
                        # Plot the final ellipsoid fully opaque
                        surface!(ax, ellipsoid.X, ellipsoid.Y, ellipsoid.Z, 
                                 color = idx_ell, colormap = (:berlin, 0.5), colorrange = (1, Nμ))
                else
                        # Plot the intermediate ellipsoids semi-transparent
                        surface!(ax, ellipsoid.X, ellipsoid.Y, ellipsoid.Z, 
                                 color = idx_ell, colormap = (:berlin, 0.01), colorrange = (1, Nμ), rasterize = 4)
                end
        end

        # Export the figure
        savefig("ellipsoid/$idx_sim.png", fig)
end

# Plot the coefficients' solutions across the parameter sweep 
function plot_mean_solutions(θm, μ_set, idx_sim, stepsize, steps)
        # Create the figure
        fig = Figure(; size = (1200, 400))

        ax1 = Axis(fig[1,1])
        ax2 = Axis(fig[1,2], title = "dt=$stepsize, Nt=$steps")
        ax3 = Axis(fig[1,3])

        # Plot the coefficients 
        lines!(ax1, μ_set, θm[:,1], color = :red, linewidth = 3.0) 
        lines!(ax2, μ_set, θm[:,2], color = :red, linewidth = 3.0)
        lines!(ax3, μ_set, θm[:,3], color = :red, linewidth = 3.0)

        # Export the figure
        savefig("solutions/$idx_sim.png", fig)
end

# Plot spectral decay
function plot_spectral_decay(λ, λc, θm, μ_set, idx_sim, stepsize, steps)
        # Create the figure
        fig = Figure(; size = (1200, 400))

        ax1L = Axis(fig[1,1])
        ax2L = Axis(fig[1,2], title = "dt=$stepsize, Nt=$steps")
        ax3L = Axis(fig[1,3])

        ax1R = Axis(fig[1,1]; yaxisposition = :right)
        hidespines!(ax1R)
        ax1R.rightspinevisible = true
        hidexdecorations!(ax1R)
        linkxaxes!(ax1L, ax1R)

        ax2R = Axis(fig[1,2]; yaxisposition = :right)
        hidespines!(ax2R)
        ax2R.rightspinevisible = true
        hidexdecorations!(ax2R)
        linkxaxes!(ax2L, ax2R)

        ax3R = Axis(fig[1,3]; yaxisposition = :right)
        hidespines!(ax3R)
        ax3R.rightspinevisible = true
        hidexdecorations!(ax3R)
        linkxaxes!(ax3L, ax3R)

        # Plot the mean of the sample
        scatter!(ax3L, μ_set, θm, markersize = 10, color = :black)

        # Plot the eigenvalues
        scatter!(ax1L, μ_set, λ[:,1], markersize = 15, color = (:red,0.5)) 
        scatter!(ax1R, μ_set, λc[:,1], markersize = 15, color = :blue) 
        scatter!(ax2L, μ_set, λ[:,2], markersize = 15, color = (:red,0.5)) 
        scatter!(ax2R, μ_set, λc[:,2], markersize = 15, color = :blue) 
        scatter!(ax3L, μ_set, λ[:,3], markersize = 15, color = (:red,0.5))
        scatter!(ax3R, μ_set, λc[:,3], markersize = 15, color = :blue) 

        # Export the figure
        savefig("spectrum/$idx_sim.png", fig)
end

# Plot the error decays with the bifurcation parameter 
function plot_error_decay(μ_set, outliers_error, interquartile_error, median_error, idx_sim, stepsize, steps)
        # Create the figure
        fig = Figure()
        ax = Axis(fig[1,1], yscale = log10, title = "dt=$stepsize, Nt=$steps")

        # Plot the error decays
        band!(ax, μ_set, outliers_error[:,1], outliers_error[:,2], color = (:steelblue, 0.25))
        lines!(ax, μ_set, outliers_error[:,1], color = :steelblue, linewidth = 3.0)
        lines!(ax, μ_set, outliers_error[:,2], color = :steelblue, linewidth = 3.0)
        lines!(ax, μ_set, interquartile_error[:,1], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax, μ_set, interquartile_error[:,2], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax, μ_set, median_error, color = :red, linewidth = 3.0)

        # Export the figure
        savefig("error_decay/$idx_sim.png", fig)
end

# Plot the error heatmap
function plot_errormap(error_matrix)
        # Create the figure
        fig = Figure()
        ax = Axis(fig[1,1], 
                  xscale = log10, yscale = log10,
                  xticks = [1e-3, 1e-2, 1e-1],
                  yticks = [1e+3, 1e+4, 1e+5]
                 )

        # Reshape the error vector into a matrix
        error = reshape(error_matrix[:,3], length(dt), length(Nt))

        # Transport the coordinate axes in logscale
        logedges(v) = (l = log10.(v); exp10.([l[1] - (l[2]-l[1])/2;
                                      (l[1:end-1] .+ l[2:end]) ./ 2;
                                      l[end] + (l[end]-l[end-1])/2]))
        x, y = logedges(dt), logedges(Nt)

        # Plot the heatmap of the error and add a colorbar
        hm = heatmap!(ax, x, y, error, colorrange = (minimum(error), maximum(error)), colorscale = log10)
        Colorbar(fig[1,2], hm)

        # Plot the location of the error minimizer in the heatmap
        idx_min = argmin(error_matrix[:,3])
        scatter!(ax, error_matrix[idx_min,1], error_matrix[idx_min,2], marker = :star5, color = :yellow, markersize = 20)

        # Export the figure
        savefig("error_matrix.png", fig)
end
