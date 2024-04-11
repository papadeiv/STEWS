import numpy as np
from matplotlib import pyplot as plt
from Process import Process

class Estimator:
    def __init__(self, process):
        self.process = process
        self.timeseries = process.realizations 
        if self.process.is_detrended:
            self.detrended = self.process.detrended

    def variance(self, width:int, is_detrended=False):
        print(self.process.x0)
        return

    def covariance(self, width:int, is_detrended=False):
        return

    def autocorrelation(self, width:int, lag:int, is_detrended=False):
        return

    def spectrum(self, width:int, is_detrended=False):
        return
