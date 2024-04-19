import numpy as np
from matplotlib import pyplot as plt

class Estimator:
    def __init__(self, timeseries):
        # Initialise the timeseries
        self.timeseries = timeseries
        # Initialise the size of the timeseries (number of observations)
        self.Nt = self.timeseries.size
        # Initialise default value of the width of the sliding window
        self.M = 0

    def mean(self, width:int):
        # Update the width of the sliding window relative to the overall length of the series
        self.M = (np.ceil(self.Nt*(width/100))).astype(int)
        # Initialise empty array
        self.smean = np.empty(self.Nt-self.M)
        # Estimate the mean across the sliding window
        for t in range(self.Nt-self.M):
            self.smean[t] = (1.0/self.M)*np.sum(self.timeseries[t:t+self.M])

        # Plot the sample mean
        self.plot('mean')
        return

    def variance(self, width:int):
        # Update the width of the sliding window relative to the overall length of the series
        self.M = (np.ceil(self.Nt*(width/100))).astype(int)
        # Compute the mean 
        self.mean(width)
        # Initialise empty array
        self.var = np.empty(self.Nt-self.M)
        # Estimate the variance across the sliding window
        for t in range(self.Nt-self.M):
            self.var[t] = (1.0/self.M)*np.sum(np.square(self.timeseries[t:t+self.M]-self.smean[t]))

        # Plot the sample variance
        self.plot('variance')
        return

    def covariance(self, width:int):
        # Update the width of the sliding window relative to the overall length of the series
        self.M = (np.ceil(self.Nt*(width/100))).astype(int)
        # Initialise empty array
        self.acvf = np.empty(self.Nt-self.M)
        # Estimate the acvf across the sliding window
        return

    def autocorrelation(self, width:int, lag=1):
        # Update the width of the sliding window relative to the overall length of the series
        self.M = (np.ceil(self.Nt*(width/100))).astype(int)
        # Compute the mean
        self.mean(width)
        # Initialise empty array
        self.acf = np.empty(self.Nt-self.M)
        # Estimate the acf across the sliding window
        for t in range(self.Nt-self.M):
            running = self.timeseries[t+lag:t+self.M]
            lagging = self.timeseries[t:(t+self.M)-lag]
            self.acf[t] = (np.sum(running - self.smean[t]))*(np.sum(lagging - self.smean[t]))/(np.sum(np.square(self.timeseries[t:t+self.M])))

        # Plot the autocorrelation coefficient
        self.plot('autocorr')
        return

    def spectrum(self, width:int):
        # Update the width of the sliding window relative to the overall length of the series
        self.M = (np.ceil(self.Nt*(width/100))).astype(int)
        # Initialise empty array
        self.spectrum = np.empty(self.Nt-self.M)
        # Estimate the spectrum across the sliding window
        return

    def plot(self, indicator):
        # Create the figure
        fig = plt.figure(figsize=[12.8,9.6], dpi=200, layout='tight')
        if indicator=='mean':
            plt.plot(np.linspace(self.M/10, (self.Nt)/10, self.smean.size), self.smean)
            plt.plot(np.linspace(self.M/10, (self.Nt)/10, self.timeseries[self.M:self.Nt].size), self.timeseries[self.M:self.Nt], color = 'black')

        elif indicator=='variance':
            ax1 = plt.subplot2grid((2,4), (0,0), colspan=4, fig=fig)
            ax1.plot(np.linspace(self.M/10, (self.Nt)/10, self.var.size), self.var)
            ax2 = plt.subplot2grid((2,4), (1,0), colspan=4, fig=fig)
            ax2.plot(np.linspace(self.M/10, (self.Nt)/10, self.timeseries[self.M:self.Nt].size), self.timeseries[self.M:self.Nt], color = 'black')

        elif indicator=='covariance':
            ax1 = plt.subplot2grid((2,4), (0,0), colspan=4, fig=fig)
            ax1.plot(np.linspace(self.M/10, (self.Nt)/10, self.acvf.size), self.acvf)
            ax2 = plt.subplot2grid((2,4), (1,0), colspan=4, fig=fig)
            ax2.plot(np.linspace(self.M/10, (self.Nt)/10, self.timeseries[self.M:self.Nt].size), self.timeseries[self.M:self.Nt], color = 'black')

        elif indicator=='autocorr':
            ax1 = plt.subplot2grid((2,4), (0,0), colspan=4, fig=fig)
            ax1.plot(np.linspace(self.M/10, (self.Nt)/10, self.acf.size), self.acf)
            ax2 = plt.subplot2grid((2,4), (1,0), colspan=4, fig=fig)
            ax2.plot(np.linspace(self.M/10, (self.Nt)/10, self.timeseries[self.M:self.Nt].size), self.timeseries[self.M:self.Nt], color = 'black')

        else:
            ax1 = plt.subplot2grid((2,4), (0,0), colspan=4, fig=fig)
            ax1.plot(np.linspace(self.M/10, (self.Nt)/10, self.smean.size), self.smean)
            ax2 = plt.subplot2grid((2,4), (1,0), colspan=4, fig=fig)
            ax2.plot(np.linspace(self.M/10, (self.Nt)/10, self.timeseries[self.M:self.Nt].size), self.timeseries[self.M:self.Nt], color = 'black')

        # Plot the statistical indicator
        plt.show()
        return
