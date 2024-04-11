#######################################################################################################################################
#      Implementing a beta process (https://en.wikipedia.org/wiki/Beta_distribution)
#######################################################################################################################################

import numpy as np
from scipy.special import beta 
from Process import Process 

# Define the gaussian distribution to draw samples from
def density(x:float):
    a = 2.0
    b = 2.0
    return (1.0/beta(a, b))*np.power(x, a-1.0)*np.power(1.0-x, b-1.0)

# Create the stochastic process object from the sampling pdf 
process = Process(density, domain=(0,1))

# Generate the timeseries of the process with Nt realizations 
Nt = 1000
process.evolve(Nt)

# Detrend the timeseries for better analysis
process.detrend(mode='diff')

# Plot the histogram of the drawn samples
process.plot()
