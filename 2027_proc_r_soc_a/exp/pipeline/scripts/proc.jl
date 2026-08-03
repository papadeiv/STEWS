"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Ground truth of the Langevin dynamical process 
U(x, μ) = 1 + μ*x - 2*(x^2) + x^4
ρ(x, μ) = exp(-U(x, μ)/D)

# Local minimum and maximum of an arbitrary cubic polynomial of coefficients
function stationary_points(coeff)
        # Compute the potential coefficients from the drift coefficients
        θ = [-coeff[1], -coeff[2]/2, -coeff[3]/3]

        # Compute the stationary points
        xs = +(1/(3*θ[3]))*(sqrt((θ[2])^2 - 3*θ[1]*θ[3]) - θ[2])
        xu = -(1/(3*θ[3]))*(sqrt((θ[2])^2 - 3*θ[1]*θ[3]) + θ[2])
        return (xu, xs)
end

# Error array 
error = Matrix{Float64}(undef, convert(Integer, Ne), 4)
