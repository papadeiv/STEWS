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
        # Solve the ensemble problem
        x0 = [get_equilibria(f, μ0, domain=[-10,10]).stable[2], μ0]
        ensemble = evolve(f, η, Λ, x0, endparameter=μf, stepsize=dt, particles=Ne)
        t = ensemble.time
        μ = ensemble.parameter
        x = ensemble.state[1]

        # -----------------------------#
        #     Statistical analysis     #
        # -----------------------------#
         
        sample_mean = Matrix{Real}(undef, convert(Integer, Ne), 3)
        sample_std = Matrix{Real}(undef, convert(Integer, Ne), 3)
        sample_skew = Matrix{Real}(undef, convert(Integer, Ne), 3)
        @showprogress for (idx_sol, solution) in enumerate(ensemble.state)
                # Compute the residuals
                x1 = detrend(solution; alg = "mean").residuals
                x2 = detrend(solution; alg = "linear", timestamps = t).residuals
                x3 = detrend(solution; alg = "emd", n_modes = 0).residuals

                # Compute the sample mean of the residuals
                sample_mean[idx_sol, 1] = mean(x1)
                sample_mean[idx_sol, 2] = mean(x2)
                sample_mean[idx_sol, 3] = mean(x3)

                # Compute the sample variance of the residuals
                sample_std[idx_sol, 1] = std(x1)
                sample_std[idx_sol, 2] = std(x2)
                sample_std[idx_sol, 3] = std(x3)

                # Compute the sample skewness of the residuals
                sample_skew[idx_sol, 1] = skewness(x1)
                sample_skew[idx_sol, 2] = skewness(x2)
                sample_skew[idx_sol, 3] = skewness(x3)
        end

        # Compute the deviations from the OUP stationary density (relative errors)
        xs = get_equilibria(f, μf, domain=[-10,10]).stable[2]
        sample_mean = abs.(sample_mean .- transpose(0.0.*ones(3)))
        sample_std = abs.(sample_std .- transpose(sqrt(σ^2/θ(xs)).*ones(3)))#./sqrt(σ^2/θ(xs))
        sample_skew = abs.(sample_skew .- transpose(0.0.*ones(3)))

        println("---------- k=1 (sample mean) ----------")
        println("   centering: ($(minimum(sample_mean[:,1])), $(mean(sample_mean[:,1])), $(maximum(sample_mean[:,1])))")
        println("  linear fit: ($(minimum(sample_mean[:,2])), $(mean(sample_mean[:,2])), $(maximum(sample_mean[:,2])))")
        println("         EMD: ($(minimum(sample_mean[:,3])), $(mean(sample_mean[:,3])), $(maximum(sample_mean[:,3])))")
        println("-------- k=2 (sample variance) --------")
        println("   centering: ($(minimum(sample_std[:,1])), $(mean(sample_std[:,1])), $(maximum(sample_std[:,1])))")
        println("  linear fit: ($(minimum(sample_std[:,2])), $(mean(sample_std[:,2])), $(maximum(sample_std[:,2])))")
        println("         EMD: ($(minimum(sample_std[:,3])), $(mean(sample_std[:,3])), $(maximum(sample_std[:,3])))")
        println("-------- k=3 (sample skewness) --------")
        println("   centering: ($(minimum(sample_skew[:,1])), $(mean(sample_skew[:,1])), $(maximum(sample_skew[:,1])))")
        println("  linear fit: ($(minimum(sample_skew[:,2])), $(mean(sample_skew[:,2])), $(maximum(sample_skew[:,2])))")
        println("         EMD: ($(minimum(sample_skew[:,3])), $(mean(sample_skew[:,3])), $(maximum(sample_skew[:,3])))")

        # Compute the optimal number of bins (Scott's rule, 1985)
        Nt = length(t)
        Nb = convert(Integer, ceil(abs(maximum(x)-minimum(x))/(3.49*std(x)*(Nt)^(-1.0/3.0))))

        # Plot the drifiting solutions
        lines!(ax1, t, x, color = :black, linewidth = 1.0)
        lines!(ax2, t, x, color = :black, linewidth = 1.0)
        lines!(ax3, t, x, color = :black, linewidth = 1.0)

        #-----------------------------#
        #     Sample mean removal     #
        #-----------------------------#
         
        # Compute and plot the trend 
        decomposition = detrend(x; alg = "mean")
        lines!(ax1, t, decomposition.trend, color = :red, linewidth = 3.0)
        # Compute and plot the distribution of the residuals 
        bins, pdf = fit_distribution(decomposition.residuals, n_bins = Nb) 
        barplot!(ax4, bins, pdf, color = (:red,0.5), strokecolor = :black, strokewidth = 1)
        # Compute and plot the stationary OUP distribution
        domain = LinRange(-0.12,0.12,1000)
        lines!(ax4, domain, [p(x,xs) for x in domain], color = :red, linewidth = 3.0)
        # Format the axes
        ax1.limits = (t[1],t[end],0.95*minimum(x),1.05*maximum(x))
        ax1.xticks = [t[1],t[end]]
        ax1.yticks = [0.95*minimum(x),1.05*maximum(x)]
        ax4.limits = (-0.12,0.12,0,14)
        ax4.xticks = [-0.12,0,0.12]
        ax4.yticks = [0,14]

        #---------------------------#
        #     Linear detrending     #
        #---------------------------#
        
        # Compute and plot the trend 
        decomposition = detrend(x; alg = "linear", timestamps = t)
        lines!(ax2, t, decomposition.trend, color = :green, linewidth = 3.0)
        # Compute and plot the distribution of the residuals 
        bins, pdf = fit_distribution(decomposition.residuals, n_bins = Nb) 
        barplot!(ax5, bins, pdf, color = (:green,0.5), strokecolor = :black, strokewidth = 1)
        # Compute and plot the stationary OUP distribution
        lines!(ax5, domain, [p(x,xs) for x in domain], color = :green, linewidth = 3.0)
        # Format the axes
        ax2.limits = (t[1],t[end],0.95*minimum(x),1.05*maximum(x))
        ax2.xticks = [t[1],t[end]]
        ax2.yticks = [0.95*minimum(x),1.05*maximum(x)]
        ax5.limits = (-0.12,0.12,0,14)
        ax5.xticks = [-0.12,0,0.12]
        ax5.yticks = [0,14]

        #--------------------------------------#
        #     Empirical mode decomposition     #
        #--------------------------------------#
         
        # Compute and plot the trend 
        decomposition = detrend(x; alg = "emd", n_modes = 1)
        lines!(ax3, t, decomposition.trend, color = :blue, linewidth = 3.0)
        # Compute and plot the distribution of the residuals 
        bins, pdf = fit_distribution(decomposition.residuals, n_bins = Nb) 
        barplot!(ax6, bins, pdf, color = (:blue,0.5), strokecolor = :black, strokewidth = 1)
        # Compute and plot the stationary OUP distribution
        lines!(ax6, domain, [p(x,xs) for x in domain], color = :blue, linewidth = 3.0)
        # Format the axes
        ax3.limits = (t[1],t[end],0.95*minimum(x),1.05*maximum(x))
        ax3.xticks = [t[1],t[end]]
        ax3.yticks = [0.95*minimum(x),1.05*maximum(x)]
        ax6.limits = (-0.12,0.12,0,14)
        ax6.xticks = [-0.12,0,0.12]
        ax6.yticks = [0,14]

        # Export the figure
        savefig("residuals_distribution.pdf", fig)
end

# Execute the main
main()
