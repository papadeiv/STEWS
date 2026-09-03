"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Define the (relative) sliding window size
window_size = 0.10::Float64

# Converts the sliding window problem into an ensemble one 
function preprocess_solution(timestamps, timeseries, width; verbose = true)
        # Assemble the sliding window
        window = build_window(length(timeseries), width)
        Nw = window.size 
        Ns = window.strides

        # Convert the sliding window subseries into an ensemble of timeseries
        if verbose
                printstyled("Converting the truncated sample path to an ensemble of ", Ns," trajectories of ", Nw, " steps\n"; bold=true, underline=true, color=:light_blue)
        end
        timesteps = [@view timestamps[n:(n+Nw-1)] for n in 1:Ns] 
        ensemble = [@view timeseries[n:(n+Nw-1)] for n in 1:Ns] 

        # Export the ensemble problem 
        return (
                timesteps = timesteps,
                trajectories = ensemble 
               ) 
end


function pullback_rate(t)
        z = -(μ0 + ε*t)/ε^(2/3)
        airy_ratio(z) = z > 0 ? airyaiprimex(z)/airyaix(z) : airyaiprime(z)/airyai(z)
        return 2*ε^(1/3)*airy_ratio(z)
end
