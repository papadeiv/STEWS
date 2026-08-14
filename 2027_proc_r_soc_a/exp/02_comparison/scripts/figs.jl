"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1200, 400))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

ax1 = Axis(fig[1,1], 
           title = "From LLS solutions",
           xlabel = "stepsize", 
           ylabel = "n. of timesteps", 
           xscale = log10, 
           yscale = log10,
           #xticks = [1e-3, 1e-2, 1e-1],
           #yticks = [1e+3, 1e+4, 1e+5]
          )
ax2 = Axis(fig[1,2], 
           title = "From AR(1) fit",
           xlabel = "stepsize", 
           xscale = log10, 
           yscale = log10,
           #xticks = [1e-3, 1e-2, 1e-1],
          )
ax3 = Axis(fig[1,3], 
           title = "From OUP variance",
           xlabel = "stepsize", 
           xscale = log10, 
           yscale = log10,
           #xticks = [1e-3, 1e-2, 1e-1],
          )
