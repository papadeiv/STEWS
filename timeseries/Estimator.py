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

    def covariance(self, width:int, lag=1):
        # Deal with the case of no-lag autocovariance
        if lag==0:
            self.variance(width)

        else:
            # Update the width of the sliding window relative to the overall length of the series
            self.M = (np.ceil(self.Nt*(width/100))).astype(int)
            # Initialise empty array
            self.acvf = np.empty(self.Nt-self.M)
            # Estimate the acvf across the sliding window
            for t in range(self.Nt-self.M):
                # Extract the current (running) sub-series and compute its mean
                running = self.timeseries[t+lag:t+self.M]
                mu_run = (1.0/(self.M - lag))*np.sum(running) 
                # Extract the past (lagging) sub-series and compute its mean
                lagging = self.timeseries[t:(t+self.M)-lag]
                mu_lag = (1.0/(self.M - lag))*np.sum(lagging) 
                # Compute the covariance
                self.acvf[t] = (1.0/(self.M - lag))*np.sum(np.multiply(running-mu_run, lagging-mu_lag))

            self.plot('covariance')

        return

    def autocorrelation(self, width:int, lag=1):
        # Update the width of the sliding window relative to the overall length of the series
        self.M = (np.ceil(self.Nt*(width/100))).astype(int)
        # Compute the acvf at lag-k (numerator)
        self.covariance(width, lag=lag)
        # Compute the acvf at lag-0 (denominator)
        self.covariance(width, lag=0)
        # Initialise empty array
        self.acf = np.empty(self.Nt-self.M)
        # Estimate the acf across the sliding window
        for t in range(self.Nt-self.M):
            self.acf[t] = self.acvf[t]/self.var[t] 

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
            plt.plot(np.linspace(self.M/10, (self.Nt)/10, self.smean.size), self.smean, label='mean')
            plt.plot(np.linspace(self.M/10, (self.Nt)/10, self.timeseries[self.M:self.Nt].size), self.timeseries[self.M:self.Nt], color = 'black')
            plt.legend()

        elif indicator=='variance':
            ax1 = plt.subplot2grid((2,4), (0,0), colspan=4, fig=fig)
            ax1.plot(np.linspace(self.M/10, (self.Nt)/10, self.var.size), self.var, label='variance')
            ax2 = plt.subplot2grid((2,4), (1,0), colspan=4, fig=fig)
            ax2.plot(np.linspace(self.M/10, (self.Nt)/10, self.timeseries[self.M:self.Nt].size), self.timeseries[self.M:self.Nt], color = 'black')
            ax1.legend()

        elif indicator=='covariance':
            ax1 = plt.subplot2grid((2,4), (0,0), colspan=4, fig=fig)
            ax1.plot(np.linspace(self.M/10, (self.Nt)/10, self.acvf.size), self.acvf, label='autocovariance')
            ax2 = plt.subplot2grid((2,4), (1,0), colspan=4, fig=fig)
            ax2.plot(np.linspace(self.M/10, (self.Nt)/10, self.timeseries[self.M:self.Nt].size), self.timeseries[self.M:self.Nt], color = 'black')
            ax1.legend()

        elif indicator=='autocorr':
            ax1 = plt.subplot2grid((2,4), (0,0), colspan=4, fig=fig)
            ax1.plot(np.linspace(self.M/10, (self.Nt)/10, self.acf.size), self.acf, label='autocorrelation')
            ax2 = plt.subplot2grid((2,4), (1,0), colspan=4, fig=fig)
            ax2.plot(np.linspace(self.M/10, (self.Nt)/10, self.timeseries[self.M:self.Nt].size), self.timeseries[self.M:self.Nt], color = 'black')
            ax1.legend()

        else:
            ax1 = plt.subplot2grid((2,4), (0,0), colspan=4, fig=fig)
            ax1.plot(np.linspace(self.M/10, (self.Nt)/10, self.smean.size), self.smean)
            ax2 = plt.subplot2grid((2,4), (1,0), colspan=4, fig=fig)
            ax2.plot(np.linspace(self.M/10, (self.Nt)/10, self.timeseries[self.M:self.Nt].size), self.timeseries[self.M:self.Nt], color = 'black')

        # Plot the statistical indicator
        plt.show()
        return
