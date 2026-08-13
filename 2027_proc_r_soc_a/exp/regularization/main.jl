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
        for (parameter_index, μ) in enumerate(μ_set[1])
                # Solve the ensemble problem 
                x0 = [sqrt(-μ), μ]
                ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

                # Loop over the regularization coefficient values
                misfit = Matrix{Float64}(undef, convert(Integer, Ne), length(α_set))
                θ_error = Matrix{Float64}(undef, convert(Integer, Ne), length(α_set))
                mean_misfit = Vector{Float64}(undef, length(α_set))
                mean_θ_error = Vector{Float64}(undef, length(α_set))
                printstyled("μ = $(μ): solving the RLLS problems\n"; bold=true, underline=true, color=:light_blue)
                @showprogress for (reg_index, α) in enumerate(α_set)
                        # Loop over the ensemble solutions
                        θ = Matrix{Float64}(undef, convert(Integer, Ne), 3)
                        for (solution_index, solution) in enumerate(ensemble.state)
                                # Check for tipping
                                if length(solution) == length(ensemble.time) 
                                        # Solve the LLS problem
                                        results = solve_lls(solution, α=α)
                                        θ[solution_index,:] = results.solution
                                        misfit[solution_index, reg_index] = results.misfit 

                                        # Compare the coefficients against the ground truth
                                        θ_error[solution_index, reg_index] = norm(results.solution - θt(μ))/norm(θt(μ))
                                else
                                        # Interrupt the execution and throw an error
                                        throw("Trajectory n. $solution_index at parameter value μ = $μ has tipped")
                                end
                        end

                        # Compute the mean misfit and the mean coefficient error of the ensemble
                        mean_misfit[reg_index] = mean(misfit[:,reg_index])
                        mean_θ_error[reg_index] = mean(θ_error[:,reg_index])
                end

                # Loop over the ensemble solutions
                for n in 1:length(ensemble.state)
                        # Plot the misfit and error dependance on the penalizer 
                        lines!(ax1, α_set, misfit[n,:], color = (:red,0.25), linewidth = 4.0)
                        lines!(ax2, α_set, θ_error[n,:], color = (:red,0.25), linewidth = 4.0)
                end

                # Plot the ensemble mean misfit and the error dependence on the penalizer
                lines!(ax1, α_set, mean_misfit, color = :blue, linewidth = 4.0)
                lines!(ax2, α_set, mean_θ_error, color = :blue, linewidth = 4.0)
        end

        # Export the figure
        savefig("marozov_discrepancy_principle.pdf", fig)
end

# Execute the main
main()
