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
           xlabel = L"\mu", 
           ylabel = "leading eigenvalue", 
           title = "From LLS solutions"
          )

ax2 = Axis(fig[1,2], 
           xlabel = L"\mu", 
           title = "From AR(1) fit"
          )

ax3 = Axis(fig[1,3], 
           xlabel = L"\mu", 
           title = "From OUP variance"
          )
