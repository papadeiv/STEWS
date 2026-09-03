"""
Numerical approximations of functions associated to the dynamics and independent from the data.

Author: Davide Papapicco 
Affil: U. of Auckland
Date: 03-08-2026
"""

function get_error(ground_truth::Function, approximation::Function, interval::Tuple, parameter::Real)
        # Define the integrands
        μ = parameter
        U(x) = ground_truth(x, μ)
        V(x) = approximation(x)
        f(x) = (U(x) - V(x))^2

        # Compute the quadrature approximation of the L2-norm of U
        integral = IntegralProblem((x, μ) -> U(x), interval)
        quadrature = solve(integral, QuadGKJL(; order=100); maxiters=1000, reltol=1e-16, abstol=1e-16)
        U2 = quadrature.u

        # Compute the quadrature approximation of the L2-norm of f
        integral = IntegralProblem((x, μ) -> f(x), interval)
        quadrature = solve(integral, QuadGKJL(; order=100); maxiters=1000, reltol=1e-16, abstol=1e-16)
        f2 = quadrature.u

        # Return the relative error
        return f2/U2 
end

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
