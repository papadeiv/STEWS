#######################################################################################################################################
#        
#       
#       
#       
#######################################################################################################################################

import numpy as np
from sklearn.linear_model import LinearRegression
from matplotlib import pyplot as plt
from Sampler import Sampler

class Process:
    def __init__(self, density, domain=(float('-inf'),float('inf'))):
        # Initialise non-uniform random variate generator (sampling) object
        self.sampler = Sampler(density, domain=domain)
        # Initialise default initial condition for the timeseries
        self.x0 = 0.0
        # Initialise default boolean variable for plotting the detrended timeseries
        self.is_detrended = False
        # Initialise default number of rows for the timeseries' plot
        self.rows = 2
        return

    def setIC(self, x0:float):
        # Specify initial condition for the timeseries
        self.x0 = x0
        return

    def evolve(self, Nt:int, drift=None, season=None):
        # Generate non-uniform variate from the sampling object according to the input pdf
        self.sampler.sample(Nt)
        # Draw uniformly at random from those variates to generate the steps of the stochastic process
        steps = self.sampler.draws
        draws = np.random.choice(Nt, Nt, replace=False)
        # Evolve the stochastic process from the initial condition through the random steps
        self.realizations = np.zeros([Nt,1])
        self.realizations[0] = self.x0
        for n in np.arange(Nt):
            if drift==None and season==None:
                self.realizations[n] = self.realizations[n-1] + steps[draws[n]]

            elif season==None:
                self.realizations[n] = drift(n/200) + steps[draws[n]]

            elif drift==None:
                self.realizations[n] = season(n/200) + steps[draws[n]] 

            else:
                self.realizations[n] = drift(n/200) + season(n) + steps[draws[n]] 

        return

    def plot(self):
        # Plot the histogram of the drawn samples
        self.sampler.plot()
        # Create a figure
        fig = plt.figure(figsize=[12.8,9.6], dpi=200, layout='tight')
        # Plot the histogram of the realizations of the stochastic process 
        ax1 = plt.subplot2grid((2,4), (0,0), colspan=1, rowspan=self.rows, fig=fig)
        ax1.hist(self.sampler.draws, bins=50, density=True, edgecolor="black", orientation='horizontal')
        # Overlap the histogram plot with the real input pdf
        x = np.linspace(np.min(self.sampler.draws), np.max(self.sampler.draws), self.realizations.size)
        ax1.plot(self.sampler.density.pdf(x), x, alpha=0.5, color = 'red', lw=1.5)
        # Plot the timeseries realizations
        ax2 = plt.subplot2grid((2,4), (0,1), colspan=3, rowspan=self.rows, fig=fig)
        ax2.plot(np.linspace(0, (self.realizations.size)/200, self.realizations.size), self.realizations, color = 'black')
        # Plot the detrended timeseries
        if self.is_detrended:
            ax3 = plt.subplot2grid((2,4), (1,1), colspan=3, rowspan=self.rows, sharex=ax2, fig=fig)
            ax3.plot(np.linspace(0, (self.realizations.size)/200, self.realizations.size), self.detrended, color = 'black')
            # Overlap the trend estimation with the original timeseries
            ax2.plot(np.linspace(0, (self.realizations.size)/200, self.realizations.size), self.trend, alpha=0.5, color = 'green')

        # Show the figure
        plt.show()

    def detrend(self, mode='EMD'):
        # Detrending by curve-fitting
        if mode=='fit':
            # Initialise discrete timesteps
            x = np.arange(self.realizations.size) 
            x = np.reshape(x, (self.realizations.size, 1))
            # Initialise the model object
            model = LinearRegression()
            # Fit the model by least-squares minimization
            model.fit(x,self.realizations)
            # Generate the trend
            self.trend = model.predict(x)
            # Detrend the original timeseries
            self.detrended = np.array(self.realizations - self.trend)

        # Detrending by first-order differencing
        elif mode=='diff':
            # Initialised detrended array
            self.detrended = np.zeros(self.realizations.shape)
            # Initialise dummy array for the trend
            self.trend = np.empty(self.realizations.shape)
            # Detrend the timeseries
            for n in np.arange(self.realizations.size-1):
                self.detrended[n] = self.realizations[n] - self.realizations[n+1]

        # Detrending by Empirical Mode Decomposition (EMD)
        else:
            print('to be implemented')

        # Update boolean variable for plotting the detrended timeseries
        self.is_detrended = True
        # Update the number of rows in the timeseries' plot
        self.rows -= 1
        return
