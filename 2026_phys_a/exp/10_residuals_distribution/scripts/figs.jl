"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1300, 650), figure_padding = (30,30,30,30))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

# Set boundaries of the potential range
upper = 4
lower = -8

# Ramped timeseries
ax1 = Axis(fig[1,1],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           ylabel = L"X_t",
           xlabelsize = labels,
           ylabelsize = labels,
           ytickformat = values -> ["$(trunc(value, digits=2))" for value in values],
           xticklabelsize = labels,
           yticklabelsize = labels,
          )
ax2 = Axis(fig[1,2],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = false,
           xticklabelsvisible = true,
           yticklabelsvisible = false,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabelsize = labels,
           ylabelsize = labels,
           ytickformat = values -> ["$(trunc(value, digits=2))" for value in values],
           xticklabelsize = labels,
           yticklabelsize = labels,
          )
ax3 = Axis(fig[1,3],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = false,
           xticklabelsvisible = true,
           yticklabelsvisible = false,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabelsize = labels,
           ylabelsize = labels,
           ytickformat = values -> ["$(trunc(value, digits=2))" for value in values],
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

# Sample mean removal
ax4 = Axis(fig[2,1],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = false,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           ylabel = L"\text{distribution}",
           xlabelsize = labels,
           ylabelsize = labels,
           xtickformat = values -> ["$(trunc(value, digits=2))" for value in values],
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

# Linear detrend
ax5 = Axis(fig[2,2],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = false,
           xticklabelsvisible = true,
           yticklabelsvisible = false,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabelsize = labels,
           ylabelsize = labels,
           xtickformat = values -> ["$(trunc(value, digits=2))" for value in values],
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

# Empirical mode decomposition
ax6 = Axis(fig[2,3],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = false,
           xticklabelsvisible = true,
           yticklabelsvisible = false,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabelsize = labels,
           ylabelsize = labels,
           xtickformat = values -> ["$(trunc(value, digits=2))" for value in values],
           xticklabelsize = labels,
           yticklabelsize = labels,
          )
