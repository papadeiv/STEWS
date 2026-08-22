"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Scalar potentials of the conservative system 
U(x, μ) = - μ*x + 2*x^2 - x^4                   # Ground truth
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)      # Arbitrary 

# Flag for tipped states
tipped = false

# Create the figure
fig = Figure(; size = (1200, 400))
ax1 = Axis(fig[1,1], 
           title = "μ=$(μ_set[1])", 
           xlabel = "stepsize", 
           ylabel = "n. of timesteps", 
           xscale = log10, 
           yscale = log10,
           xticks = [1e-2, 1e-1],
           yticks = [1e+3, 1e+4]
          )
ax2 = Axis(fig[1,2], 
           title = "μ=$(μ_set[2])", 
           xlabel = "stepsize", 
           xscale = log10, 
           yscale = log10,
           yticklabelsvisible = false,
           xticks = [1e-2, 1e-1],
           yticks = [1e+3, 1e+4]
          )
ax3 = Axis(fig[1,3], 
           title = "μ=$(μ_set[3])", 
           xlabel = "stepsize", 
           xscale = log10, 
           yscale = log10,
           yticklabelsvisible = false,
           xticks = [1e-2, 1e-1],
           yticks = [1e+3, 1e+4]
          )
axes = [ax1, ax2, ax3]

# Solve the LLS problem
function solve_lls(solution, dt)
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

# Perform the sensitivity analysis of the generated samples
function analysis()
        # Parameter sweep loop
        printstyled("Analyzing the samples\n"; bold=true, underline=true, color=:light_blue)
        err_min, err_max = 10, 0
        relative_error = Matrix{Real}(undef, Ns, 3)
        @showprogress for (idx_μ, μ) in enumerate(μ_set)
                # Compute the ground truth of the modified escape rate 
                xs = (get_equilibria(f, μ, domain=[-10,10])).stable[2]
                xu = (get_equilibria(f, μ, domain=[-10,10])).unstable[1]
                ΔU = abs(U(xu, μ) - U(xs, μ))
                escape = exp(-ΔU)

                # Simulation parameters loop
                for (idx_sim, (stepsize, steps)) in enumerate(Iterators.product(dt,Nt))
                        # Import the data the ensemble median of the modified escape rate
                        ews = median(readin("$idx_μ/$stepsize/$steps.csv"))

                        # Compute the relative error between the median estimate and the ground truth
                        relative_error[idx_sim, idx_μ] = abs(escape-ews)/abs(escape)
                end

                # Extract the minimum and maximum relative errors
                err_min > minimum(relative_error[:,idx_μ]) && (err_min = minimum(relative_error[:,idx_μ]))
                err_max < maximum(relative_error[:,idx_μ]) && (err_max = maximum(relative_error[:,idx_μ]))
        end

        # Parameter sweep loop
        local hm
        @showprogress for (idx_μ, μ) in enumerate(μ_set)
                # Transformation of the coordinate axes in logscale so that it does not interpolate negative values
                logedges(v) = (l = log10.(v); exp10.([l[1] - (l[2]-l[1])/2;
                              (l[1:end-1] .+ l[2:end]) ./ 2;
                              l[end] + (l[end]-l[end-1])/2]))
                x, y = logedges(dt), logedges(Nt)

                # Plot the heatmap of the relative errors and the minimizers
                error_map = reshape(relative_error[:,idx_μ], length(dt), length(Nt))
                hm = heatmap!(axes[idx_μ], x, y, error_map, colorrange = (err_min, err_max), colormap = Reverse(:Paired_11))
        end

        # Add the colorbar 
        Colorbar(fig[1,4],
                 size = 25,
                 label = "relative error",
                 hm,
                )

        # Export the figure
        savefig("13_sensitivity_analysis.pdf", fig)
end
