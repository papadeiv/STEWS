"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Estimate the return rate EWS across a sliding window
function analyse(time, state, width, idx)
        # Assemble the ensemble of windowed subseries
        window = build_window(time, state, width)

        # Loop over the window's strides
        data = Matrix{Float64}(undef, length(window.timesteps), 4)
        for (idx_x, x) in enumerate(window.trajectories)
                # Append timestamp to the return rate
                data[idx_x,1] = (window.timesteps[idx_x])[end]

                # Estimate the return rate in the current window
                data[idx_x,2] = (estimate_return_rate(x, dt, alg = "LLS"))^2
                data[idx_x,3] = (estimate_return_rate(x, dt, alg = "AC1"))^2
                data[idx_x,4] = (estimate_return_rate(x, dt, alg = "OUP", σ=σ))^2
        end

        # Export the data 
        writeout(data, "extrapolation/$idx.csv")

        return length(window.timesteps)
end

function extrapolate(time, strides)
        # Loop over the sample paths
        tipping_times = Matrix{Float64}(undef, convert(Integer, Ne), 3)
        return_rate_1 = Matrix{Float64}(undef, strides, convert(Integer, Ne) + 1)
        return_rate_2 = Matrix{Float64}(undef, strides, convert(Integer, Ne) + 1)
        return_rate_3 = Matrix{Float64}(undef, strides, convert(Integer, Ne) + 1)
        for idx_sol in 1:convert(Integer, Ne)
                # Import the data 
                data = readin("extrapolation/$idx_sol.csv")

                # Append estimated return rate in ensemble data structure
                return_rate_1[:,1] = data[:,1]
                return_rate_2[:,1] = data[:,1]
                return_rate_3[:,1] = data[:,1]
                return_rate_1[:,idx_sol+1] = data[:,2]
                return_rate_2[:,idx_sol+1] = data[:,3]
                return_rate_3[:,idx_sol+1] = data[:,4]

                # Linear extrapolation of the return rate to find the tipping time
                coeff = estimate_tipping_time(data[:,1], data[:,2]) # LLS
                tipping_times[idx_sol,1] = -coeff[1]/coeff[2]
                coeff = estimate_tipping_time(data[:,1], data[:,3]) # AC1 
                tipping_times[idx_sol,2] = -coeff[1]/coeff[2]
                coeff = estimate_tipping_time(data[:,1], data[:,4]) # OUP 
                tipping_times[idx_sol,3] = -coeff[1]/coeff[2]
        end

        # Plot the tipping time distribution
        #hist!(ax2, tipping_times[:,1], bins=10, color = (:orange,0.5), normalization = :pdf)

        # Initialize data structures for the analysis
        lower_envelope = Matrix{Real}(undef, strides, 3)
        first_quartile = Matrix{Real}(undef, strides, 3)
        median_estimate = Matrix{Real}(undef, strides, 3)
        third_quartile = Matrix{Real}(undef, strides, 3)
        upper_envelope = Matrix{Real}(undef, strides, 3)

        # Loop over the strides of the estimated return rate
        for idx in 1:size(return_rate_1,1)
                # Compute the lower envelope of the three estimates 
                lower_envelope[idx,1] = (extrema(return_rate_1[idx,2:end]))[1]
                lower_envelope[idx,2] = (extrema(return_rate_2[idx,2:end]))[1]
                lower_envelope[idx,3] = (extrema(return_rate_3[idx,2:end]))[1]

                # Compute the 1st quartile of the three estimates
                first_quartile[idx,1] = quantile(return_rate_1[idx,2:end], 0.25)
                first_quartile[idx,2] = quantile(return_rate_2[idx,2:end], 0.25)
                first_quartile[idx,3] = quantile(return_rate_3[idx,2:end], 0.25)

                # Compute the median (2nd quartile) of the three estimates
                median_estimate[idx,1] = median(return_rate_1[idx,2:end])
                median_estimate[idx,2] = median(return_rate_2[idx,2:end])
                median_estimate[idx,3] = median(return_rate_3[idx,2:end])

                # Compute the 3rd quartile of the three estimates
                third_quartile[idx,1] = quantile(return_rate_1[idx,2:end], 0.75)
                third_quartile[idx,2] = quantile(return_rate_2[idx,2:end], 0.75)
                third_quartile[idx,3] = quantile(return_rate_3[idx,2:end], 0.75)

                # Compute the upper envelope of the three estimates
                upper_envelope[idx,1] = (extrema(return_rate_1[idx,2:end]))[2]
                upper_envelope[idx,2] = (extrema(return_rate_2[idx,2:end]))[2]
                upper_envelope[idx,3] = (extrema(return_rate_3[idx,2:end]))[2]
        end

        # Plot the range 
        band!(ax2, return_rate_1[:,1], lower_envelope[:,1], upper_envelope[:,1], color = (:steelblue, 0.25))
        band!(ax3, return_rate_2[:,1], lower_envelope[:,2], upper_envelope[:,2], color = (:steelblue, 0.25))
        band!(ax4, return_rate_3[:,1], lower_envelope[:,3], upper_envelope[:,3], color = (:steelblue, 0.25))

        # Plot the envelopes
        lines!(ax2, return_rate_1[:,1], lower_envelope[:,1], color = :steelblue, linewidth = 2.0)
        lines!(ax2, return_rate_1[:,1], upper_envelope[:,1], color = :steelblue, linewidth = 2.0)
        lines!(ax3, return_rate_2[:,1], lower_envelope[:,2], color = :steelblue, linewidth = 2.0)
        lines!(ax3, return_rate_2[:,1], upper_envelope[:,2], color = :steelblue, linewidth = 2.0)
        lines!(ax4, return_rate_3[:,1], lower_envelope[:,3], color = :steelblue, linewidth = 2.0)
        lines!(ax4, return_rate_3[:,1], upper_envelope[:,3], color = :steelblue, linewidth = 2.0)

        # Plot the interquartile range
        band!(ax2, return_rate_1[:,1], first_quartile[:,1], third_quartile[:,1], color = :palegreen)
        band!(ax3, return_rate_2[:,1], first_quartile[:,2], third_quartile[:,2], color = :palegreen)
        band!(ax4, return_rate_3[:,1], first_quartile[:,3], third_quartile[:,3], color = :palegreen)

        # Plot the first and third quantiles 
        lines!(ax2, return_rate_1[:,1], first_quartile[:,1], color = :darkgreen, linewidth = 1.0)
        lines!(ax2, return_rate_1[:,1], third_quartile[:,1], color = :darkgreen, linewidth = 1.0)
        lines!(ax3, return_rate_2[:,1], first_quartile[:,2], color = :darkgreen, linewidth = 1.0)
        lines!(ax3, return_rate_2[:,1], third_quartile[:,2], color = :darkgreen, linewidth = 1.0)
        lines!(ax4, return_rate_3[:,1], first_quartile[:,3], color = :darkgreen, linewidth = 1.0)
        lines!(ax4, return_rate_3[:,1], third_quartile[:,3], color = :darkgreen, linewidth = 1.0)

        # Plot the median estimates
        lines!(ax2, return_rate_1[:,1], median_estimate[:,1], color = :brown1, linewidth = 2.0)
        lines!(ax3, return_rate_2[:,1], median_estimate[:,2], color = :brown1, linewidth = 2.0)
        lines!(ax4, return_rate_3[:,1], median_estimate[:,3], color = :brown1, linewidth = 2.0)

        # Plot the linear extrapolation of the median return rate
        coeff = estimate_tipping_time(return_rate_1[:,1], median_estimate[:,1]) # LLS
        median_tipping_time = -coeff[1]/coeff[2]
        timerange = LinRange(return_rate_1[end,1], median_tipping_time, 100)
        lines!(ax2, timerange, coeff[2].*timerange .+ coeff[1], color = :brown1, linewidth = 3.0)
        hlines!(ax2, 0, color = (:black, 0.5), linewidth = 2.0)
        scatter!(ax2, median_tipping_time, 0, markersize = 20, color = :brown1, strokewidth = 2.0)

        coeff = estimate_tipping_time(return_rate_2[:,1], median_estimate[:,2]) # AC1 
        median_tipping_time = -coeff[1]/coeff[2]
        timerange = LinRange(return_rate_2[end,1], median_tipping_time, 100)
        lines!(ax3, timerange, coeff[2].*timerange .+ coeff[1], color = :brown1, linewidth = 3.0)
        hlines!(ax3, 0, color = (:black, 0.5), linewidth = 2.0)
        scatter!(ax3, median_tipping_time, 0, markersize = 20, color = :brown1, strokewidth = 2.0)

        coeff = estimate_tipping_time(return_rate_3[:,1], median_estimate[:,3]) # LLS
        median_tipping_time = -coeff[1]/coeff[2]
        timerange = LinRange(return_rate_3[end,1], median_tipping_time, 100)
        lines!(ax4, timerange, coeff[2].*timerange .+ coeff[1], color = :brown1, linewidth = 3.0)
        hlines!(ax4, 0, color = (:black, 0.5), linewidth = 2.0)
        scatter!(ax4, median_tipping_time, 0, markersize = 20, color = :brown1, strokewidth = 2.0)
end
