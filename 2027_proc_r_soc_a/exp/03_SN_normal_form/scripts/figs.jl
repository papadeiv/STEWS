"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1000, 500))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

ax1 = Axis(fig[1,1], 
           xlabel = L"\mu", 
           ylabel = L"(Df(x^*))^2", 
           title = "From LLS solutions"
          )

ax2 = Axis(fig[1,2], 
           xlabel = L"\mu", 
           title = "From KBR solutions"
          )

ax3 = Axis(fig[1,3], 
           xlabel = L"\mu", 
           title = "From Bayesian estimate"
          )
