"""
Approximation methods to reconstruct the scalar potential of a conservative system from the empirical distribution of a sample path.

Author: Davide Papapicco Affil: U. of Auckland
Date: 26-09-2025
"""

# Approximate the normalization constant of Boltzmann-like pdf using a generic potential
function normalise(f::Function, parameter::Real; accuracy=1e-16)
        # Define the integrand
        μ = parameter
        ρ(x) = f(x, μ)

        # Solve the definite integral by using adaptive Gauss-Kronrod quadrature
        integral = IntegralProblem((x, μ) -> ρ(x), (-Inf, Inf))
        quadrature = solve(integral, QuadGKJL(; order=100); maxiters=1000, reltol=accuracy, abstol=accuracy)

        # Compute the normalisation constant
        N = 1.0::Float64/(quadrature.u)
        return N
end

# Approximate the normalization constant of a Boltzmann-like pdf using a cubic potential
function normalise(f::Function, interval::Tuple; accuracy=1e-8, order=100, maxiters=1000)
        # Use order=20000 and maxiters=10000 for the NLLS method
         
        # Solve the definite integral by using adaptive Gauss-Kronrod quadrature
        integral = IntegralProblem((x, _) -> f(x), interval)
        quadrature = solve(integral, QuadGKJL(; order=order); maxiters=maxiters, reltol=accuracy, abstol=accuracy)

        # Compute the normalisation constant
        N = 1.0::Float64/(quadrature.u)
        return N
end

# Numerical approximation of the relative error of the difference potential (U-V)(x) in the L2-norm
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
