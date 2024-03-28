#######################################################################################################################################
#      Implementing a simple random walker under gaussian distribution (https://en.wikipedia.org/wiki/Random_walk#Gaussian_random_walk)
#######################################################################################################################################

import numpy as np
from Sampler import Sampler

# Define the gaussian distribution to draw samples from
def density(x:float):
    mu = 0.0
    sigma = 1.00
    return (1.0/sigma*np.sqrt(2*np.pi))*(np.exp(-0.5*np.square((x-mu)/sigma)))

# Create non-uniform sampling object
sampler = Sampler(density, domain=(-10,10))
# Draw samples from the input pdf
sampler.sample(10000)
# Plot the histogram of the drawn samples
sampler.plot()
