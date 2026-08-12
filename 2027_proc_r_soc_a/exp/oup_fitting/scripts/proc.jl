"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Flag for tipped states
tipped = false

# Aribtrary quartic potential
V(x, θ) = θ[1]*x + θ[2]*(x^2)

# Ground truth of the Langevin dynamics and its Taylor approximation
U(x, μ) = 1 + μ*x - 2*(x^2) + x^4
Ux(x, μ) = μ -4*x + 4*x^3
U2x(x, μ) = -4 + 12*x^2
U3x(x, μ) = + 24*x
T(x, x0, μ) = U(x0, μ) + Ux(x0, μ)*(x - x0) + (U2x(x0, μ)/2)*(x - x0)^2 + (U3x(x0, μ)/6)*(x - x0)^3

# Taylor approximation in monomial basis
c0(x0, μ) = U(x0, μ) - Ux(x0, μ)*x0 + (U2x(x0, μ)/2)*(x0^2) - (U3x(x0, μ)/6)*(x0^3)
c1(x0, μ) = Ux(x0, μ) - U2x(x0, μ)*x0 + (U3x(x0, μ)/2)*(x0^2)
c2(x0, μ) = U2x(x0, μ)/2 - (U3x(x0, μ)/2)*x0
c3(x0, μ) = U3x(x0, μ)/6
Tm(x, x0, μ) = c0(x0, μ) + c1(x0, μ)*x + c2(x0, μ)*(x^2) + c3(x0, μ)*(x^3)

# Initialize data structures for the analysis
median_solutions = [zeros(Nμ) for _ in 1:2]
interquartile_solutions = [zeros(Nμ, 2) for _ in 1:2]
outliers_solutions = [zeros(Nμ, 2) for _ in 1:2]
spectrum = Matrix{Real}(undef, Nμ, 2)

# Perform the statistical analysis of the generated samples
function analysis()
        # Loop over the parameter values
        printstyled("Analyzing the samples\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for (idx_μ, μ) in enumerate(μ_set)
                # Import the data and convert it from drift to potential coefficients
                data = readin("solutions/$idx_μ.csv")
                θ = hcat(-data[:,1], -data[:,2]./2)

                for (idx_col, col) in enumerate(eachcol(θ))
                        median_solutions[idx_col][idx_μ] = median(col)
                        interquartile_solutions[idx_col][idx_μ,:] = quantile(col, [0.25, 0.75])
                        outliers_solutions[idx_col][idx_μ,:] = collect(extrema(col))
                end
        end

        # Plot and export the solutions, error and spectral decay figures
        plot_solutions(μ_set, outliers_solutions, interquartile_solutions, median_solutions)
end
