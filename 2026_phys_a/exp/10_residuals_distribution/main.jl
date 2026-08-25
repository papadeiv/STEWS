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
        ensemble = evolve(f, η, Λ, x0, endparameter=μf, stepsize=dt)
        t = ensemble.time
        μ = ensemble.parameter
        x = ensemble.state[1]

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
        xs = get_equilibria(f, μf, domain=[-10,10]).stable[2]
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
