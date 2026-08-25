"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1200, 800))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

# Set boundaries of the potential range
upper = 4
lower = -8

# Top row (simulation parameters)
ax1 = Axis(fig[1,1], 
           title = "μ=$(μ_set[1])", 
           xscale = log10, 
           yscale = log10,
           xticks = [1e-2, 1e-1],
           #yticks = [1e+3, 1e+4]
          )
ax2 = Axis(fig[1,2], 
           title = "μ=$(μ_set[2])", 
           xscale = log10, 
           yscale = log10,
           yticklabelsvisible = false,
           xticks = [1e-2, 1e-1],
           #yticks = [1e+3, 1e+4]
          )
ax3 = Axis(fig[1,3], 
           title = "μ=$(μ_set[3])", 
           xscale = log10, 
           yscale = log10,
           yticklabelsvisible = false,
           xticks = [1e-2, 1e-1],
           #yticks = [1e+3, 1e+4]
          )
top_axes = [ax1, ax2, ax3]

# Bottom row (system parameters)
ax4 = Axis(fig[2,1], 
           xscale = log10, 
           yscale = log10,
           xticks = [1e-3, 1e-1],
           #yticks = [1e-3, 1e-5]
          )
ax5 = Axis(fig[2,2], 
           xscale = log10, 
           yscale = log10,
           yticklabelsvisible = false,
           xticks = [1e-3, 1e-1],
           #yticks = [1e-3, 1e-5]
          )
ax6 = Axis(fig[2,3], 
           xscale = log10, 
           yscale = log10,
           yticklabelsvisible = false,
           xticks = [1e-3, 1e-1],
           #yticks = [1e-3, 1e-5]
          )
bottom_axes = [ax4, ax5, ax6]
