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

ax1 = Axis(fig[1,1],
           title = L"\mu(t) = -1 + \varepsilon\cdot t",
           titlesize = labels,
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = true,
           xticklabelsvisible = false,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           #xticks = [μ0, μf],
           #yticks = [-1,0,1],
           xtickwidth = border,
           ytickwidth = border,
           ylabel = L"\text{state}",
           ylabelsize = labels,
           yticklabelsize = labels,
          )

ax2 = Axis(fig[1,2],
           title = L"\mu(t) = -0.5",
           titlesize = labels,
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = false,
           xticklabelsvisible = false,
           yticklabelsvisible = false,
           xtickalign = 1,
           ytickalign = 1,
           #xticks = [μ0, μf],
           #yticks = [-1,0,1],
           xtickwidth = border,
           ytickwidth = border,
          )

ax3 = Axis(fig[1,3],
           title = L"\mu(t) = -0.1",
           titlesize = labels,
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = false,
           xticklabelsvisible = false,
           yticklabelsvisible = false,
           xtickalign = 1,
           ytickalign = 1,
           #xticks = [μ0, μf],
           #yticks = [-1,0,1],
           xtickwidth = border,
           ytickwidth = border,
          )

ax4 = Axis(fig[2,1],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           #xticks = [μ0, μf],
           #yticks = [-1,0,1],
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\text{time}",
           ylabel = L"\text{ews}",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

ax5 = Axis(fig[2,2],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = false,
           xticklabelsvisible = true,
           yticklabelsvisible = false,
           xtickalign = 1,
           ytickalign = 1,
           #xticks = [μ0, μf],
           #yticks = [-1,0,1],
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\text{time}",
           xlabelsize = labels,
           xticklabelsize = labels,
          )

ax6 = Axis(fig[2,3],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = false,
           xticklabelsvisible = true,
           yticklabelsvisible = false,
           xtickalign = 1,
           ytickalign = 1,
           #xticks = [μ0, μf],
           #yticks = [-1,0,1],
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\text{time}",
           xlabelsize = labels,
           xticklabelsize = labels,
          )

top_axes = [ax1, ax2, ax3]
bottom_axes = [ax4, ax5, ax6]
