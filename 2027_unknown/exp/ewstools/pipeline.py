import ewstools
import numpy as np
import matplotlib.pyplot as plt

from ewstools.models import simulate_ricker

# Generate the sample path
series = simulate_ricker(tmax=500, F=[0,2.7])

# Obtain the quasi-stationary residuals
ts = ewstools.TimeSeries(data=series, transition=440)
ts.detrend(method='Lowess', span=0.2)

# Compute the EWS of the residuals
ts.compute_var(rolling_window=0.5)
ts.compute_auto(lag=1, rolling_window=0.5)
ts.compute_auto(lag=2, rolling_window=0.5)
ts.compute_ktau()

"""
print(vars(ts))
print(ts.state)
print(ts.transition)
print(ts.ews)
print(ts.ktau)
"""

# Plot the timeseries
series.plot();
plt.show()
