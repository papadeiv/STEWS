"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/figs.jl")

# Define the main algorithm
function main()
        # Loop over the parameter values
        printstyled("Generating the samples\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for (idx_μ, μ) in enumerate(μ_set)
                # Solve the ensemble problem 
                x0 = [sqrt(-μ), μ]
                ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

                # Loop over the regularization coefficient values
                estimates = Matrix{Float64}(undef, convert(Integer, Ne), 3)
                for (idx_sol, solution) in enumerate(ensemble.state)
                        # Check for tipping
                        if length(solution) == length(ensemble.time) 
                                # Solve the LLS problem
                                u = solution .- x0[1]
                                estimates[idx_sol,1] = ((estimate_parameters(u, dt))[2]).^2

                                # Lamouroux et al. 2009
                                estimates[idx_sol,2] = (estimate_parameters_2009(solution, dt).α).^2

                                # Willers et al. 2024
                                estimates[idx_sol,3] = (estimate_parameters_2024(solution, dt).α).^2
                        else
                                # Interrupt the execution and throw an error
                                throw("Trajectory n. $solution_index at parameter value μ = $μ has tipped")
                        end
                end

                # Export the data
                writeout(estimates, "solutions/$idx_μ.csv")
        end

        # Perform the statistical analysis for the current simulation setup
        analysis()
end

# Execute the main
main()
