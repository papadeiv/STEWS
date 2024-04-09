#######################################################################################################################################
#       This class implements a non-uniform random variate generator (https://en.wikipedia.org/wiki/Inverse_transform_sampling) 
#       to draw samples from a specific probability density function (pdf) provided as input by the user. The sampling method
#       that has been used is the inverse transformation method (https://en.wikipedia.org/wiki/Inverse_transform_sampling) which
#       is implemented by the stats module in SciPy and doesn't require the knowledge of the cumulative distribution (CDF)
#######################################################################################################################################

from Sampler import Sampler

class Process:
    def __init__(self):
        return

    def setIC(self, x0:float):
        return

    def evolve(self, Nt:int):
        return

    def detrend(self, mode='EMD'):
        return
