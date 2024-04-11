#######################################################################################################################################
#      Implementing a simple random walker under gaussian distribution (https://en.wikipedia.org/wiki/Random_walk#Gaussian_random_walk)
#######################################################################################################################################

import numpy as np
from Process import Process 
from Estimator import Estimator

# Define the gaussian distribution to draw samples from
def density(x:float):
    mu = 0.0
    sigma = 1.00
    return (1.0/(sigma*np.sqrt(2*np.pi)))*(np.exp(-0.5*np.square((x-mu)/sigma)))

# Create the stochastic process object from the sampling pdf 
random_walker = Process(density, domain=(-10,10))

# Generate the timeseries of the process with Nt realizations 
Nt = 1000
random_walker.evolve(Nt)

# Detrend the timeseries for better analysis
random_walker.detrend(mode='diff')

# Plot the histogram of the drawn samples
random_walker.plot()

# Estimate statistical indicators of the timeseries
estimator = Estimator(random_walker)
estimator.variance(5)
