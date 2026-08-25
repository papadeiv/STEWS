"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

using Distributions
using StatsBase
using HypothesisTests
using Random
using Statistics

# Ornstein-Uhlenbeck process' (OUP) stationary density
θ(x) = 12*x^2 - 4
p(x,xs) = sqrt(θ(xs)/(2*π*σ^2))*exp(-(θ(xs)*x^2)/(2*σ^2))
