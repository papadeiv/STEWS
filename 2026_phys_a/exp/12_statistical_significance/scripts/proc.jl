"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Scalar potential of the conservative system 
U(x, μ) = + μ*x + (1.0/3.0)*x^3                # Ground truth
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)      # Reconstructed 

# Define the (relative) sliding window size
window_size = 0.10::Float64

# Converts the sliding window problem into an ensemble one 
function preprocess_solution(timestamps, parameters, timeseries, width)
        # Assemble the sliding window
        window = build_window(length(timeseries), width)
        Nw = window.size
        Ns = window.strides

        # Convert the sliding window subseries into an ensemble of timeseries
        timesteps = [@view timestamps[n:(n+Nw-1)] for n in 1:Ns] 
        ramps = [@view parameters[n:(n+Nw-1)] for n in 1:Ns] 
        ensemble = [@view timeseries[n:(n+Nw-1)] for n in 1:Ns] 

        # Export the ensemble problem 
        return (
                timesteps = timesteps,
                parameters = ramps,
                trajectories = ensemble 
               ) 
end

# Solve the LLS problem
function solve_lls(solution)
        # Define the observation vectors 
        Xn = solution[1:end-1]
        Y  = (solution[2:end] .- solution[1:end-1])./dt

        # Assemble the model matrix
        A = hcat(ones(length(Xn)), Xn, Xn.^2)

        # Solve the (linear) least-squares problem
        β = A\Y

        # Compute the coefficients of the potential
        θ = [-β[1], -β[2]/2, -β[3]/3]
        return θ 
end
 
# Compute the modified escape EWS
function compute_ews(θ)
        # Compute estimated stable and unstable equilibria of the cubic approximation
        xs = +(1/(3*θ[3]))*(sqrt((θ[2])^2 - 3*θ[1]*θ[3]) - θ[2])
        xu = -(1/(3*θ[3]))*(sqrt((θ[2])^2 - 3*θ[1]*θ[3]) + θ[2])

        # Compute the modified escape EWS 
        ΔV = abs(V(xu, θ) - V(xs, θ))
        escape = exp(-ΔV)

        # Return the EWS 
        return escape 
end

# Perform the statistical analysis of the generated samples
function analysis(data_idx)
        # Import the data 
        data = readin("ews/$data_idx.csv")
        t = (readin("ews/ground_truth.csv"))[:,1]

        # Initialize data structures for the analysis
        lower_envelope = Vector{Real}(undef, length(t))
        first_quartile = Vector{Real}(undef, length(t))
        median_estimate = Vector{Real}(undef, length(t))
        third_quartile = Vector{Real}(undef, length(t))
        upper_envelope = Vector{Real}(undef, length(t))

        # Loop over the timesteps
        for (row_idx, row) in enumerate(eachrow(data))
                # Compute the statistics of the data
                lower_envelope[row_idx] = (extrema(row))[1]
                first_quartile[row_idx] = quantile(row, 0.25)
                median_estimate[row_idx] = median(row)
                third_quartile[row_idx] = quantile(row, 0.75)
                upper_envelope[row_idx] = (extrema(row))[2]
        end

        # Plot the range 
        band!(bottom_axes[data_idx], t, lower_envelope, upper_envelope, color = (:steelblue, 0.25))

        # Plot the envelopes
        lines!(bottom_axes[data_idx], t, lower_envelope, color = :steelblue, linewidth = 2.0)
        lines!(bottom_axes[data_idx], t, upper_envelope, color = :steelblue, linewidth = 2.0)

        # Plot the interquartile range
        band!(bottom_axes[data_idx], t, first_quartile, third_quartile, color = :palegreen)

        # Plot the first and third quantiles 
        lines!(bottom_axes[data_idx], t, first_quartile, color = :darkgreen, linewidth = 1.0)
        lines!(bottom_axes[data_idx], t, third_quartile, color = :darkgreen, linewidth = 1.0)

        # Plot the median estimates
        lines!(bottom_axes[data_idx], t, median_estimate, color = :brown1, linewidth = 2.0)
end
