# Plot the coefficients' solutions across the parameter sweep 
function plot_solutions(μ_set, outliers_solutions, interquartile_solutions, median_solutions, outliers_variance, interquartile_variance, median_variance)
        # Create the figure
        fig = Figure(; size = (1200,600))
        ax1 = Axis(fig[1,1], xlabel = L"\mu", ylabel = L"\tilde{c}_0")
        ax2 = Axis(fig[1,2], xlabel = L"\mu", ylabel = L"\tilde{c}_1")
        ax3 = Axis(fig[1,3], xlabel = L"\mu", ylabel = L"\tilde{c}_2")
        ax4 = Axis(fig[1,4], xlabel = L"\mu", ylabel = L"\tilde{c}_3")
        ax5 = Axis(fig[2,1], xlabel = L"\mu", ylabel = L"U(x_0)")
        ax6 = Axis(fig[2,2], xlabel = L"\mu", ylabel = L"U^{\,'}(x_0)")
        ax7 = Axis(fig[2,3], xlabel = L"\mu", ylabel = L"U^{\,''}(x_0)/2")
        ax8 = Axis(fig[2,4], xlabel = L"\mu", ylabel = L"U^{\,'''}(x_0)/6")

        # Plot the coefficient solutions at different values of the bifucation paramater
        # c1
        band!(ax6, μ_set, outliers_solutions[1][:,1], outliers_solutions[1][:,2], color = (:steelblue, 0.25))
        lines!(ax6, μ_set, outliers_solutions[1][:,1], color = :steelblue, linewidth = 3.0)
        lines!(ax6, μ_set, outliers_solutions[1][:,2], color = :steelblue, linewidth = 3.0)
        lines!(ax6, μ_set, interquartile_solutions[1][:,1], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax6, μ_set, interquartile_solutions[1][:,2], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax6, μ_set, median_solutions[1], color = :red, linewidth = 3.0)
        # c2
        band!(ax7, μ_set, outliers_solutions[2][:,1], outliers_solutions[2][:,2], color = (:steelblue, 0.25))
        lines!(ax7, μ_set, outliers_solutions[2][:,1], color = :steelblue, linewidth = 3.0)
        lines!(ax7, μ_set, outliers_solutions[2][:,2], color = :steelblue, linewidth = 3.0)
        lines!(ax7, μ_set, interquartile_solutions[2][:,1], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax7, μ_set, interquartile_solutions[2][:,2], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax7, μ_set, median_solutions[2], color = :red, linewidth = 3.0)
        # c3
        band!(ax8, μ_set, outliers_solutions[3][:,1], outliers_solutions[3][:,2], color = (:steelblue, 0.25))
        lines!(ax8, μ_set, outliers_solutions[3][:,1], color = :steelblue, linewidth = 3.0)
        lines!(ax8, μ_set, outliers_solutions[3][:,2], color = :steelblue, linewidth = 3.0)
        lines!(ax8, μ_set, interquartile_solutions[3][:,1], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax8, μ_set, interquartile_solutions[3][:,2], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax8, μ_set, median_solutions[3], color = :red, linewidth = 3.0)

        # Loop over the parameter values
        new_μ_set = LinRange(-1.2, 1.5396, Nμ)
        taylor_coeff = Matrix{Real}(undef, Nμ, 4)
        monomial_coeff = Matrix{Real}(undef, Nμ, 4)
        for (idx_μ, μ) in enumerate(new_μ_set)
                # Compute the stable equilibrium at the current parameter value
                x0 = (get_equilibria(f, μ, domain=[-10,10])).stable[2]

                # Compute the coefficients of the Taylor expansion of U
                taylor_coeff[idx_μ,1] = U(x0, μ)
                taylor_coeff[idx_μ,2] = 0.0 
                taylor_coeff[idx_μ,3] = U2x(x0, μ)/2
                taylor_coeff[idx_μ,4] = U3x(x0, μ)/6

                # Compute the monomimial coefficients of the Taylor expansion of U
                monomial_coeff[idx_μ,1] = c0(x0, μ)
                monomial_coeff[idx_μ,2] = c1(x0, μ)
                monomial_coeff[idx_μ,3] = c2(x0, μ)
                monomial_coeff[idx_μ,4] = c3(x0, μ)

                #=
                # Plot the potential and its Taylor approximations
                fig2 = Figure(; size = (1200, 600))
                ax21 = Axis(fig2[1,1], limits = (-2, 2, -2, 2), title = "μ=$(μ_set[idx_μ])")
                ax22 = Axis(fig2[1,2], limits = (-2, 2, -2, 2), title = "μ=$(μ_set[idx_μ])")
                domain = LinRange(-2,2,1000)
                lines!(ax21, domain, [U(x, μ) for x in domain], color = :black, linewidth = 5.0)
                lines!(ax21, domain, [T(x, x0, μ) for x in domain], color = (:steelblue,0.75), linewidth = 5.0)
                lines!(ax21, domain, [Tm(x, x0, μ) for x in domain], color = :blue, linewidth = 2.0)

                band!(ax22, domain, [V(x, [outliers_solutions[1][idx_μ,1], outliers_solutions[2][idx_μ,1], outliers_solutions[3][idx_μ,1]]) for x in domain], [V(x, [outliers_solutions[1][idx_μ,2], outliers_solutions[2][idx_μ,2], outliers_solutions[3][idx_μ,2]]) for x in domain], color = (:steelblue, 0.25))
                lines!(ax22, domain, [V(x, [outliers_solutions[1][idx_μ,1], outliers_solutions[2][idx_μ,1], outliers_solutions[3][idx_μ,1]]) for x in domain], color = :steelblue, linewidth = 3.0)
                lines!(ax22, domain, [V(x, [outliers_solutions[1][idx_μ,2], outliers_solutions[2][idx_μ,2], outliers_solutions[3][idx_μ,2]]) for x in domain], color = :steelblue, linewidth = 3.0)
                lines!(ax22, domain, [V(x, [interquartile_solutions[1][idx_μ,1], interquartile_solutions[2][idx_μ,1], interquartile_solutions[3][idx_μ,1]]) for x in domain], color = :steelblue, linestyle = :dash, linewidth = 3.0)
                lines!(ax22, domain, [V(x, [interquartile_solutions[1][idx_μ,2], interquartile_solutions[2][idx_μ,2], interquartile_solutions[3][idx_μ,2]]) for x in domain], color = :steelblue, linestyle = :dash, linewidth = 3.0)
                lines!(ax22, domain, [V(x,[median_solutions[1][idx_μ], median_solutions[2][idx_μ], median_solutions[3][idx_μ]]) for x in domain], color = :red, linewidth = 3.0)

                # Export the figure
                savefig("potentials/$idx_μ.png", fig2)
                =#
        end

        # Plot the Taylor's expansion coefficients at different values of the bifurcation parameter
        lines!(ax5, new_μ_set, taylor_coeff[:,1], color = :black, linewidth = 2.0)
        lines!(ax6, new_μ_set, taylor_coeff[:,2], color = :black, linewidth = 2.0)
        lines!(ax7, new_μ_set, taylor_coeff[:,3], color = :black, linewidth = 2.0)
        lines!(ax8, new_μ_set, taylor_coeff[:,4], color = :black, linewidth = 2.0)

        # Plot the Taylor's expansion monomial coefficients at different values of the bifurcation parameter
        lines!(ax1, new_μ_set, monomial_coeff[:,1], color = :black, linewidth = 2.0)
        lines!(ax2, new_μ_set, monomial_coeff[:,2], color = :black, linewidth = 2.0)
        lines!(ax3, new_μ_set, monomial_coeff[:,3], color = :black, linewidth = 2.0)
        lines!(ax4, new_μ_set, monomial_coeff[:,4], color = :black, linewidth = 2.0)

        # Export the figure
        savefig("solutions.png", fig)

        # Create the figure
        fig = Figure(; size = (600,600))
        ax1 = Axis(fig[1,1], xlabel = L"\mu", ylabel = "inv. leading eigenvalue", title = "From LLS solutions", limits = (new_μ_set[1],new_μ_set[end],0,6))
        ax2 = Axis(fig[1,2], xlabel = L"\mu", title = "From OUP variance", limits = (new_μ_set[1],new_μ_set[end],0,6))
        ax3 = Axis(fig[2,1], xlabel = L"\mu", ylabel = L"\text{Var}[X_t]", title = "W.r.t. stationary variance", limits = (new_μ_set[1],new_μ_set[end],0.0001,0.0005))
        ax4 = Axis(fig[2,2], xlabel = L"\mu", title = "W.r.t. OUP variance", limits = (new_μ_set[1],new_μ_set[end],0.0001,0.0005))

        # Plot the leading eigenvalue reconstruction from LLS 
        band!(ax1, μ_set, outliers_solutions[2][:,1], outliers_solutions[2][:,2], color = (:steelblue, 0.25))
        lines!(ax1, μ_set, outliers_solutions[2][:,1], color = :steelblue, linewidth = 3.0)
        lines!(ax1, μ_set, outliers_solutions[2][:,2], color = :steelblue, linewidth = 3.0)
        lines!(ax1, μ_set, interquartile_solutions[2][:,1], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax1, μ_set, interquartile_solutions[2][:,2], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax1, μ_set, median_solutions[2], color = :red, linewidth = 3.0)

        # Plot the leading eigenvalue reconstruction from OUP 
        band!(ax2, μ_set, (σ^2/4).*inv.(outliers_variance[:,1]), (σ^2/4).*inv.(outliers_variance[:,2]), color = (:steelblue, 0.25))
        lines!(ax2, μ_set, (σ^2/4).*inv.(outliers_variance[:,1]), color = :steelblue, linewidth = 3.0)
        lines!(ax2, μ_set, (σ^2/4).*inv.(outliers_variance[:,2]), color = :steelblue, linewidth = 3.0)
        lines!(ax2, μ_set, (σ^2/4).*inv.(interquartile_variance[:,1]), color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax2, μ_set, (σ^2/4).*inv.(interquartile_variance[:,2]), color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax2, μ_set, (σ^2/4).*inv.(median_variance), color = :red, linewidth = 3.0)
        
        # Plot the ground truth of the Taylor expansion
        lines!(ax1, new_μ_set, taylor_coeff[:,3], color = :black, linewidth = 2.0)
        lines!(ax2, new_μ_set, taylor_coeff[:,3], color = :black, linewidth = 2.0)

        # Plot the sample variance 
        band!(ax3, μ_set, outliers_variance[:,1], outliers_variance[:,2], color = (:steelblue, 0.25))
        lines!(ax3, μ_set, outliers_variance[:,1], color = :steelblue, linewidth = 3.0)
        lines!(ax3, μ_set, outliers_variance[:,2], color = :steelblue, linewidth = 3.0)
        lines!(ax3, μ_set, interquartile_variance[:,1], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax3, μ_set, interquartile_variance[:,2], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax3, μ_set, median_variance, color = :red, linewidth = 3.0)

        # Plot the OUP variance 
        band!(ax4, μ_set, outliers_variance[:,1], outliers_variance[:,2], color = (:steelblue, 0.25))
        lines!(ax4, μ_set, outliers_variance[:,1], color = :steelblue, linewidth = 3.0)
        lines!(ax4, μ_set, outliers_variance[:,2], color = :steelblue, linewidth = 3.0)
        lines!(ax4, μ_set, interquartile_variance[:,1], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax4, μ_set, interquartile_variance[:,2], color = :steelblue, linestyle = :dash, linewidth = 3.0)
        lines!(ax4, μ_set, median_variance, color = :red, linewidth = 3.0)

        # Plot the ground truth of the OUP variance
        lines!(ax4, new_μ_set[1:end-1], ((σ^2/4).*inv.(taylor_coeff[:,3]))[1:end-1], color = :black, linewidth = 2.0)

        # Export the figure
        savefig("new_EWS?.png", fig)
end
