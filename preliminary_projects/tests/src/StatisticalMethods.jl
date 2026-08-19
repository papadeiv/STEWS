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
export normalise, get_error
export build_window, detrend, find_tipping
export fit_potential, fit_density, fit_ellipsoid
export estimate_distribution, estimate_parameters

end # module
