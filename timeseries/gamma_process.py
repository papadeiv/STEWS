#######################################################################################################################################
#      Implementing a gamma process (https://en.wikipedia.org/wiki/Gamma_process)
#######################################################################################################################################

import numpy as np
from scipy.special import gamma
from Process import Process 

# Define the gaussian distribution to draw samples from
def density(x:float):
    k = 7.5
    theta = 1.0
    return (1.0/(gamma(k)*(np.power(theta, k))))*np.power(x, k-1.0)*np.exp(-x/theta)

# Create the stochastic process object from the sampling pdf 
process = Process(density, domain=(0,20))

# Generate the timeseries of the process with Nt realizations 
Nt = 1000
process.evolve(Nt)

# Detrend the timeseries for better analysis
process.detrend(mode='fit')

# Plot the histogram of the drawn samples
process.plot()
