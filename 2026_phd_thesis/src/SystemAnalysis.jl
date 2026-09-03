module SystemAnalysis

# Import packages
using LinearAlgebra, StochasticDiffEq
using NonlinearSolve, Roots, ForwardDiff
using DocStringExtensions

# Import utility functions
include("../utils/evolution.jl")
include("../utils/equilibria.jl")
include("../utils/approximation.jl")

# Export namespaces
export evolve, get_equilibria
export get_error, fit_potential

end # module
