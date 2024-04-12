import numpy as np
from matplotlib import pyplot as plt

class Estimator:
    def __init__(self, timeseries):
        self.timeseries = timeseries 
        self.Nt = self.timeseries.size

    def variance(self, width:int):
        return

    def covariance(self, width:int):
        return

    def autocorrelation(self, width:int, lag:int):
        return

    def spectrum(self, width:int):
        return
