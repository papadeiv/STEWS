import numpy as np
from matplotlib import pyplot as plt

class Estimator:
    def __init__(self, timeseries):
        self.timeseries = timeseries 
        self.Nt = self.timeseries.size

    def mean(self, width:int):
        # Compute the size of the sliding window relative to the overall length of the series
        M = np.ceil(self.Nt*(width/100))
        M = M.astype(int)
        # Initialise empty array
        self.smean = np.empty(self.Nt-M)
        # Estimate the mean across the sliding window
        for t in range(self.Nt-M):
            self.smean[t] = (1.0/M)*np.sum(self.timeseries[t:t+M])

        # Plot the sample mean
        self.plot('mean')
        return

    def variance(self, width:int):
        # Compute the size of the sliding window relative to the overall length of the series
        M = np.ceil(self.Nt*(width/100))
        M = M.astype(int)
        # Compute the mean 
        self.mean(width)
        # Initialise empty array
        self.var = np.empty(self.Nt-M)
        # Estimate the variance across the sliding window
        for t in range(self.Nt-M):
            self.var[t] = (1.0/M)*np.sum(np.square(self.timeseries[t:t+M]-self.smean[t]))

        # Plot the sample variance
        self.plot('variance')
        return

    def covariance(self, width:int):
        # Compute the size of the sliding window relative to the overall length of the series
        M = np.ceil(self.Nt*(width/100))
        M = M.astype(int)
        # Initialise empty array
        self.acvf = np.empty(self.Nt-M)
        # Estimate the acvf across the sliding window
        return

    def autocorrelation(self, width:int, lag=1):
        # Compute the size of the sliding window relative to the overall length of the series
        M = np.ceil(self.Nt*(width/100))
        M = M.astype(int)
        # Compute the mean
        self.mean(width)
        # Initialise empty array
        self.acf = np.empty(self.Nt-M)
        # Estimate the acf across the sliding window
        for t in range(self.Nt-M):
            running = self.timeseries[t+lag:t+M]
            lagging = self.timeseries[t:(t+M)-lag]
            self.acf[t] = (np.sum(running - self.smean[t]))*(np.sum(lagging - self.smean[t]))/(np.sum(np.square(self.timeseries[t:t+M])))

        # Plot the autocorrelation coefficient
        self.plot('autocorr')
        return

    def spectrum(self, width:int):
        # Compute the size of the sliding window relative to the overall length of the series
        M = np.ceil(self.Nt*(width/100))
        M = M.astype(int)
        # Initialise empty array
        self.spectrum = np.empty(self.Nt-M)
        # Estimate the spectrum across the sliding window
        return

    def plot(self, indicator):
        # Create the figure
        plt.figure(figsize=[12.8,9.6], dpi=200, layout='tight')
        # Initialise array for the statistical indicator to be plotted
        if indicator=='mean':
            plt.plot(np.arange(self.smean.size), self.smean)

        elif indicator=='variance':
            plt.plot(np.arange(self.var.size), self.var)

        elif indicator=='covariance':
            plt.plot(np.arange(self.acvf.size), self.acvf)

        elif indicator=='autocorr':
            plt.plot(np.arange(self.acf.size), self.acf)

        else:
            plt.plot(np.arange(self.mean.size), self.mean)

        # Plot the statistical indicator
        plt.show()
        return
