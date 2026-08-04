"""
Local approximation of the probability density p and the potential V of a Langevin dynamical process. The potential V is shifted to be centered at the target local minimo x0 of U.

Author: Davide Papapicco
Affil: University of Auckland
Date: 03-08-2026
"""

# Quantiles of χ² with 3 degrees of freedom: P(‖x‖²_Σ ≤ q) = conf for a 3-D
# Gaussian.  The ellipsoid is {x : (x-μ)ᵀΣ⁻¹(x-μ) ≤ d²} with d = √q.
const χ2 = Dict(0.50   => 2.3660,
                0.6827 => 3.5267,   # 3-D analogue of "one std away from the mean"
                0.90   => 6.2514,
                0.95   => 7.8147,
                0.99   => 11.3449)

# Builds the covariance (uncertainty) ellipsoid of the estimated coefficients of the cubic potential
function fit_ellipsoid(sample; confidence = 0.95, nθ = 100, nφ = 100)
        # Compute the mean vector of the solutions
        μ = vec(mean(sample; dims = 1))

        # Compute the covariance matrix Σ and its symmetric factor L s.t. Σ = L*L^T 
        Σ = cov(sample)
        Λ = eigen(Symmetric(Matrix(Σ)))
        λ = Λ.values
        v = Λ.vectors
        L = v*Diagonal(sqrt.(max.(λ, 0)))

        # Compute the Mahalanobis distance given the confidence
        d = sqrt(χ2[confidence])

        # Parametrize the domain of the ellipsoid
        θ = LinRange(0, π, nθ)
        φ = LinRange(0, 2π, nφ)

        # Build the ellipsoid as a surface
        X = Matrix{Float64}(undef, nθ, nφ); 
        Y = similar(X); 
        Z = similar(X)
        for i in 1:nθ, j in 1:nφ
                u = [sin(θ[i])*cos(φ[j]), sin(θ[i])*sin(φ[j]), cos(θ[i])]
                p = μ .+ d.*(L*u)
                X[i, j], Y[i, j], Z[i, j] = p
        end

        return (
                mean = μ,
                covariance = Σ,
                distance = d,
                ellipsoid = (
                             X = X, 
                             Y = Y, 
                             Z = Z
                            )
               )
end

# Returns a Boltzmann-like equilibrium distribution, based on a cubic potential V, bounded on one end of the domain
function fit_density(V, D)
        # Define the Boltzmann-like functional form
        p(x) = exp(-V(x)/D)

        # Compute critical points of V numerically
        dV(x) = ForwardDiff.derivative(V, x)
        d2V(x) = ForwardDiff.derivative(dV, x)
        points = find_zeros(dV, -5, 5)
        xu = points[1]
        xs = points[2]
        if d2V(points[1]) > 0
                xs = points[1]
                xu = points[2]
        else
                xs = points[2]
                xu = points[1]
        end

        # Establish the integration interval
        interval = (-Inf, Inf)
        if xu > xs 
                interval = (-Inf, xu)
        else
                interval = (xu, +Inf)
        end

        # Compute the normalized pdf 
        N = normalise(p, interval)
        ρ(x) = N*p(x)

        return ρ, interval 
end

# Coefficients are the (regularized) linear least-squares solution of the Euler-Maruyama quasi-likelihood estimation 
function fit_potential(coeff, x0, U, μ)
        # Compute the potential coefficients from the drift coefficients
        θ = [-coeff[1], -coeff[2]/2, -coeff[3]/3]

        # Compute the stable equilibrium (center of the shift)
        xs = +(1/(3*θ[3]))*(sqrt((θ[2])^2 - 3*θ[1]*θ[3]) - θ[2])

        # Compute the shifts
        δx = x0 - xs 
        δy = U(x0, μ) - (Polynomial([0.0; θ]))(xs)

        # Define the shifted potential
        V(x) = δy + θ[1]*(x - δx) + θ[2]*(x - δx)^2 + θ[3]*(x - δx)^3

        return V
end

#=
function fit_potential(bins, distribution, degree, noise::Float64; N = nothing)
        # Compute the diffusion coefficient
        D = (noise^2)/2

        # Filter out the 0-valued entries in the distribution
        idx = findall(x -> x > 0.0, distribution)
        y = [distribution[n] for n in idx]
        x = [bins[n] for n in idx]

        # Define the normalisation constant based on user input
        if N == nothing
                # Assumption of a stationary OUP
                N=1/sqrt(2*pi*D)
        end

        # Compute the distribution on the potential using the stationarity assumption
        U = -D.*log.(y./N)

        # Solve the linear least-square problem to fit U(x) to V
        V = Polynomials.fit(x, U, degree)

        return (
                points = x,
                potential = U,
                fit = V.coeffs
               )
end

function fit_potential(timeseries; n_bins=nothing, noise=nothing, initial_guess=nothing, transformation = [1.0, 0.0, 0.0], optimiser=1e-2, attempts=1000, verbose = false)
        # Number of bins for the histogram
        if n_bins == nothing
                # Scott's rule (1985)
                n_bins = convert(Int64, ceil(abs(maximum(timeseries)-minimum(timeseries))/(3.49*std(timeseries)*(length(timeseries))^(-1.0/3.0))))   
        end
        Nb = n_bins

        # Additive noise in the SDE
        if noise == nothing
                noise = std(timeseries)
        end
        σ = noise

        # Noise in the optimiser's steps
        β = optimiser

        # Fit an empirical distribution to the timeseries data
        bins, hist = fit_distribution(timeseries, n_bins=Nb+1)

        # Find the median of the empirical distribution (for centering the polynomial weighting)
        dx = bins[2] - bins[1]
        cdf = cumsum(hist)*dx
        median_idx = findfirst(>=(0.5), cdf)

        # Initial guess for the non-linear 0-problem 
        if initial_guess == nothing
                initial_guess = (fit_potential(bins, hist, 3, σ)).fit[2:(4)]
        end

        # Define the stochastic diffusion
        D = (σ^2)/2.0::Float64

        # Compute a shift for the potential {c0} that sets V(xs)=0 to avoid numerical cancellation
        xs(μ) = (1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) - μ[2])
        c0(μ) = - μ[1]*xs(μ) - μ[2]*(xs(μ))^2 - μ[3]*(xs(μ))^3

        # Define an arbitrary cubic with the the above constraint on {c0}
        V(x, μ) = c0(μ) + μ[1]*x + μ[2]*(x^2) + μ[3]*(x^3)
        # Define the unnormalised pdf as an exponential of the abritrary cubic
        f(x, μ) = exp(-(1.0::Float64/D)*(V(x, μ)))

        # Define the normalisation constant as a function of the 4 parameters
        N(μ) = normalise(f, μ)

        # Define the target of the optimisation problem: normalised pdf
        p(x, μ) = N(μ)*exp.(-(1.0::Float64/D).*(c0(μ) .+ μ[1]*x .+ μ[2]*(x.^2) .+ μ[3]*(x.^3)))

        # Define the lower and upper bounds for the coefficients
        lower = [-45.0, -25.0, -40.0] 
        upper = [45.0, 25.0, 40.0]

        # Define the polynomial weighting
        α = transformation[1]
        β = transformation[2]
        d = convert(Int64, transformation[3])
        weights = α .+ β.*(bins .- bins[median_idx]).^d

        # First attempt to solve the non-linear least-squares problem
        try
                solution = curve_fit(p, bins, hist, weights, initial_guess, lower=lower, upper=upper).param
                return (
                        points = bins,
                        potential = nothing,
                        fit = solution
                       )
        catch e
                if isa(e, ArgumentError) && (
                                             occursin("matrix contains Infs or NaNs", e.msg) ||
                                             occursin("Initial guess must be within bounds", e.msg)
                                            )
                        # No print
                elseif isa(e, DomainError) || isa(e, LoadError)
                        # No print
                elseif isa(e, TaskFailedException) || isa(e, LinearAlgebra.SingularException)
                        # No print
                else
                        rethrow(e)
                end
        end

        # Perturbed attempts
        tries = 0
        while tries < attempts
                try
                        # Perturb the initial guess
                        perturbed_guess = initial_guess + β.*randn(3)
                        # Attempt to solve the nonlinear problem
                        solution = curve_fit(p, bins, hist, weights, perturbed_guess, lower=lower, upper=upper).param
                        return (
                                points = bins,
                                potential = nothing,
                                fit = solution
                               )
                catch e
                        if isa(e, ArgumentError) && (
                                                     occursin("matrix contains Infs or NaNs", e.msg) ||
                                                     occursin("Initial guess must be within bounds", e.msg)
                                                    )
                                tries += 1
                        elseif isa(e, DomainError) || isa(e, LoadError)
                                tries += 1
                        elseif isa(e, TaskFailedException) || isa(e, LinearAlgebra.SingularException)
                                tries += 1
                        else
                                rethrow(e)
                        end
                end
        end

        # Return the linear solution
        if verbose
                println("Curve fitting failed after $(tries) attempts.")
        end

        return (
                points = bins,
                potential = nothing,
                fit = initial_guess
               )
end
=#
