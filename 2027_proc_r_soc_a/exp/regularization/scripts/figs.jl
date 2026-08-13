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

# Marozov's discrepancy principle
ax1 = Axis(fig[1,1], limits = (nothing, nothing, nothing, nothing),
          spinewidth = border,
          xgridvisible = false,
          ygridvisible = false,
          xlabelvisible = true,
          ylabelvisible = true,
          xticklabelsvisible = true,
          yticklabelsvisible = true,
          xtickalign = 1,
          ytickalign = 1,
          xtickwidth = border,
          ytickwidth = border,
          xlabel = L"\alpha",
          ylabel = L"||\mathbf{Y} - \mathbf{A}\theta||/||\mathbf{Y}||_2",
          xlabelsize = labels,
          ylabelsize = labels,
          xticklabelsize = labels,
          yticklabelsize = labels,
         )

# Marozov's discrepancy principle
ax2 = Axis(fig[1,2], limits = (nothing, nothing, nothing, nothing),
          spinewidth = border,
          xgridvisible = false,
          ygridvisible = false,
          xtickalign = 1,
          ytickalign = 1,
          xtickwidth = border,
          ytickwidth = border,
          xlabel = L"\alpha",
          ylabel = L"||\theta-\theta_t||_2/||\theta_t||_2",
          xlabelsize = labels,
          ylabelsize = labels,
          xticklabelsize = labels,
          yticklabelsize = labels,
         )
