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

# Main algorithm 
function main()
        #----------------------------------------------#
        #     Ramped problem with positive tipping     #
        #----------------------------------------------#
        
        # Solve the ensemble problem 
        x0 = [sqrt(-μ0), μ0]
        ensemble = evolve(f, η, Λ, x0, endparameter=μf, stepsize=dt, particles=Ne)

        # Determine earliest tipping in the ensemble
        tip_idx = length(ensemble.time)
        for (sol_idx, solution) in enumerate(ensemble.state)
                # Identify the tipping and update the index accordingly
                tipping = find_tipping(solution; width = 300, threshold = 2.75, verbose=false) - 100::Integer
                tipping < tip_idx && (tip_idx = tipping)
        end

        # Truncate the timesteps and the parameter ramps
        μ = ensemble.parameter[1:tip_idx]
        t = ensemble.time[1:tip_idx]
        Nt = length(t)
        display(Nt)

        # Export the parameter ramp and the earliest tipping
        writeout(hcat(t[end], μ[end]), "solutions/tipping.csv")
        writeout(hcat(ensemble.time, ensemble.parameter), "solutions/parameters.csv")

        # Compute the lengths of windowed timeseries and strides
        window = build_window(Nt, window_size)
        Ns = window.strides 

        # Estimate the modified escape rate from the ensemble
        ground_truth = Matrix{Real}(undef, Ns, 3)
        ews = Matrix{Real}(undef, Ns, length(ensemble.state))
        printstyled("Generating the samples\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for (sol_idx, solution) in enumerate(ensemble.state)
                # Export the solution
                sol_idx ≤ 10 && (writeout(hcat(ensemble.time[1:length(solution)], solution), "solutions/1/$sol_idx.csv"))

                # Compute the residuals and extract the windowed subseries 
                residuals = detrend(solution[1:tip_idx]; alg = "emd", n_modes = 1).residuals
                subseries = preprocess_solution(t, μ, residuals, window_size)

                # Estimate the modified escate rate across the sliding windows
                for (win_idx, xw) in enumerate(subseries.trajectories)
                        # Store the timestep and parameter value at the window's end
                        sol_idx == 1 && (ground_truth[win_idx,1] = (subseries.timesteps[win_idx])[end])
                        sol_idx == 1 && (ground_truth[win_idx,2] = (subseries.parameters[win_idx])[end])

                        # Compute the ground truth
                        μw = ground_truth[win_idx,2]
                        a = sqrt(-μw)
                        b = -sqrt(-μw)
                        ΔU = abs(U(b, μw) - U(a, μw))
                        ground_truth[win_idx,3] = exp(-ΔU)

                        # Solve the LLS problem and compute the ews
                        θ = solve_lls(xw)
                        ews[win_idx,sol_idx] = compute_ews(θ)
                end
        end

        # Export the EWS timeseries
        writeout(ground_truth, "ews/ground_truth.csv")
        writeout(ews, "ews/1.csv")

        #-------------------------------------------------------------#
        #     Stationary problem with negative tipping (μ = -0.6)     #
        #-------------------------------------------------------------#
         
        # Redefine the dynamics for the stationary problem
        Λn(t) = 0.0 
        x0 = [sqrt(-μ1), μ1]

        # Solve the ensemble problem
        ensemble = evolve(f, η, Λn, x0, stepsize=dt, steps=(Nt-1), particles=Ne)
        μ = ensemble.parameter
        t = ensemble.time

        # Ensure absence of tipping 
        for (sol_idx, solution) in enumerate(ensemble.state)
                if length(solution) < length(t)
                        @goto tipped
                end
        end

        # Estimate the modified escape rate from the ensemble
        ews = Matrix{Real}(undef, Ns, length(ensemble.state))
        @showprogress for (sol_idx, solution) in enumerate(ensemble.state)
                # Export the solution
                sol_idx ≤ 10 && (writeout(hcat(ensemble.time[1:length(solution)], solution), "solutions/2/$sol_idx.csv"))

                # Compute the residuals and extract the windowed subseries 
                residuals = detrend(solution; alg = "emd", n_modes = 1).residuals
                subseries = preprocess_solution(t, μ, residuals, window_size)

                # Estimate the modified escate rate across the sliding windows
                for (win_idx, xw) in enumerate(subseries.trajectories[1:Ns])
                        θ = solve_lls(xw)
                        ews[win_idx,sol_idx] = compute_ews(θ)
                end
        end

        # Export the EWS timeseries
        writeout(ews, "ews/2.csv")

        #-------------------------------------------------------------#
        #     Stationary problem with negative tipping (μ = -0.3)     #
        #-------------------------------------------------------------#
         
        # Redefine the dynamics for the stationary problem
        x0 = [sqrt(-μ2), μ2]

        # Solve the ensemble problem
        ensemble = evolve(f, η, Λn, x0, stepsize=dt, steps=(Nt-1), particles=Ne)
        μ = ensemble.parameter
        t = ensemble.time

        # Ensure absence of tipping 
        for (sol_idx, solution) in enumerate(ensemble.state)
                if length(solution) < length(t)
                        @goto tipped
                end
        end

        # Estimate the modified escape rate from the ensemble
        ews = Matrix{Real}(undef, Ns, length(ensemble.state))
        @showprogress for (sol_idx, solution) in enumerate(ensemble.state)
                # Export the solution
                sol_idx ≤ 10 && (writeout(hcat(ensemble.time[1:length(solution)], solution), "solutions/3/$sol_idx.csv"))

                # Compute the residuals and extract the windowed subseries 
                residuals = detrend(solution; alg = "emd", n_modes = 1).residuals
                subseries = preprocess_solution(t, μ, residuals, window_size)

                # Estimate the modified escate rate across the sliding windows
                for (win_idx, xw) in enumerate(subseries.trajectories[1:Ns])
                        θ = solve_lls(xw)
                        ews[win_idx,sol_idx] = compute_ews(θ)
                end
        end

        # Export the EWS timeseries
        writeout(ews, "ews/3.csv")

        #------------------------------#
        #     Statistical analysis     #
        #------------------------------#
 
        # Loop over the three different experiments
        printstyled("Analyzing the samples\n"; bold=true, underline=true, color=:light_blue)
        for anl_idx in 1:3
                # Loop over the solutions of the ensemble problem
                for sol_idx in 1:convert(Integer, Ne)
                        if sol_idx < 10
                                # Import the data 
                                data = readin("solutions/$anl_idx/$sol_idx.csv")
                                t = data[:,1]
                                x = data[:,2]

                                # Plot the timeseries
                                lines!(top_axes[anl_idx], t, x, color = (:black,0.125), linewidth = 1.0)
                        end
                end

                # Format axes of the figure
                top_axes[anl_idx].limits = (0,20000,-1.1,1.1)
                bottom_axes[anl_idx].limits = (0,20000,0,1)

                # Perform and plot the statistical analysis on the ews
                analysis(anl_idx)
        end

        # Plot the ground truths for the top axes
        pars = readin("solutions/parameters.csv")
        tip = readin("solutions/tipping.csv")
        lines!(ax1, pars[:,1], [sqrt(-μ) for μ in pars[:,2]], color = :blue, linewidth = 2.5)
        lines!(ax1, pars[:,1], [-sqrt(-μ) for μ in pars[:,2]], color = :red, linewidth = 2.5)
        vlines!(ax1, tip[1], color = :black, linewidth = 2.5, linestyle = :dash)
        hlines!(ax2, sqrt(-μ1), color = :blue, linewidth = 2.5)
        hlines!(ax2, -sqrt(-μ1), color = :red, linewidth = 2.5)
        vlines!(ax2, tip[1], color = :black, linewidth = 2.5, linestyle = :dash)
        hlines!(ax3, sqrt(-μ2), color = :blue, linewidth = 2.5)
        hlines!(ax3, -sqrt(-μ2), color = :red, linewidth = 2.5)
        vlines!(ax3, tip[1], color = :black, linewidth = 2.5, linestyle = :dash)

        # Plot the ground truths for the bottom axes
        ews = readin("ews/ground_truth.csv")
        lines!(ax4, ews[:,1], ews[:,3], color = :black, linewidth = 2.5)
        hlines!(ax5, exp(-abs(U(-sqrt(-μ1),μ1) - U(sqrt(-μ1),μ1))), color = :black, linewidth = 2.5)
        hlines!(ax6, exp(-abs(U(-sqrt(-μ2),μ2) - U(sqrt(-μ2),μ2))), color = :black, linewidth = 2.5)

        # Export the figure
        savefig("14_false_negatives.pdf", fig)

        @label tipped
end

# Execute the main 
main()
