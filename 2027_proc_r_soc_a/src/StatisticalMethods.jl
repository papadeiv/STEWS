module StatisticalMethods

# Import packages
using LinearAlgebra, StatsBase, LsqFit
using Polynomials, Integrals
using Roots, ForwardDiff
using ProgressMeter, DocStringExtensions

# Import Python package for empirical mode decomposition
using PyCall
PyEMD = pyimport("PyEMD")

# Import utility functions
include("../utils/preprocess.jl")
include("../utils/estimation.jl")
include("../utils/approximation.jl")
include("../utils/reconstruction.jl")

using Tables, CSV, MAT, DataFrames
include("../utils/output.jl")

# Export namespaces
export build_window, detrend, find_tipping
export fit_distribution, normalise, get_error
export estimate_parameters, fit_potential, fit_density

end # module
