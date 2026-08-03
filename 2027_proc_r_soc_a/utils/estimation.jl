"""
???

Author: Davide Papapicco
Affil: U. of Auckland
Date: 30-09-2025
"""

# Approximates a probability distribution in L1 using histograms (piecewise constant functions)
function estimate_distribution(u; interval = nothing, n_bins = 200::Int64)
        # Get the range of values of the distribution
        if interval == nothing
                global u_min = u[argmin(u)]
                global u_max = u[argmax(u)]
                global range = u_max - u_min
        else
                global u_min = interval[1]
                global u_max = interval[end]
                global range = interval[end] - interval[1] 
        end

        # Define the edges of the bins of the histogram
        x = LinRange(u_min - range*0.05, u_max + range*0.05, n_bins)

        # Derive the center-points for the locations of the bins in the plot
        bins = [(x[n+1]+x[n])/2 for n in 1:(length(x)-1)]

        # Fit the histogram through the defined bins
        hist = StatsBase.fit(Histogram, u, x)

        # Normalise the histrogram to get an empirical pdf
        pdf = (LinearAlgebra.normalize(hist, mode = :pdf)).weights
                
        # Export the data
        return bins, pdf, LinearAlgebra.norm(hist)
end

# Implementation of the (regularized) linear-least squares solution of the Euler-Mauryama quasi-likelihood estimation
function estimate_parameters(timeseries, dt; α = 0.0)
        # Define the observation and data vectors 
        Xn = timeseries[1:end-1]
        Y  = (timeseries[2:end] .- timeseries[1:end-1])./dt

        # Assemble the model matrix
        A = hcat(ones(length(Xn)), Xn, Xn.^2)

        # Solve the regularized (linear) least-squares problem
        θ = (A'*A + α.*I(size(A,2)))\(A'*Y)
        return θ 
end
