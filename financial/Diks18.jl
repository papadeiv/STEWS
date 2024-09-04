using CairoMakie, Makie.Colors
using CSV, DataFrames

CairoMakie.activate!(; px_per_unit = 2)

df = DataFrame(CSV.File("./S&P500FinancialCrisis.csv"; delim=','))
time = df[!,1]
SST = df[!,2]

fig = Figure(; size = (1200, 800), backgroundcolor = :transparent)
#=
ax = Axis(fig[1,1:4],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    #limits = ((0,10),(18,27)),
    # Title
    title = "2008-09 global financial crisis",
    titlevisible = true,
    titlesize = 50,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = "days",
    xlabelvisible = true,
    xlabelsize = 50,
    xlabelcolor = :black,
    xlabelpadding = -60.0,
    xticks = [495,369,161,33],#, ["02/2008", "08/2008", "06/2009", "12/2009"]),
    xticksvisible = true,
    xticksize = 20,
    xticklabelsvisible = true,
    xticklabelsize = 50,
    #xtickformat = "{:.0f}",
    xscale = log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = "index"
    ylabelvisible = true,
    ylabelsize = 50,
    ylabelcolor = :black,
    ylabelpadding = -60.0,
    #yticks = [700,800,1300,1400],
    yticksvisible = true,
    yticksize = 20,
    yticklabelsvisible = true,
    yticklabelsize = 50,
    ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left
)
=#
ax = Axis(fig[1,1:4],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    #limits = ((0,10),(18,27)),
    # Title
    title = "2008-09 global financial crisis",
    titlevisible = true,
    titlesize = 50,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = "days",
    xlabelvisible = true,
    xlabelsize = 50,
    xlabelcolor = :black,
    xlabelpadding = -60.0,
    xticksvisible = true,
    xticksize = 20,
    xticklabelsvisible = true,
    xticklabelsize = 50,
    #xtickformat = "{:.0f}",
    xscale = identity,
    xaxisposition = :bottom,
    # y-axis
    ylabel = "index (kUSD)",
    ylabelvisible = true,
    ylabelsize = 50,
    ylabelcolor = :black,
    ylabelpadding = -60.0,
    yticks = [7,14],
    yticksvisible = true,
    yticksize = 20,
    yticklabelsvisible = true,
    yticklabelsize = 50,
    ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
    xticks = ([495,369,139,20],["02/2008", "08/2008", "07/2009", "01/2010"])
)
ax.xreversed = true
#text!(475, 7.5, text = L"\text{(f)}", align = [:right, :top], color = :black, fontsize = 70)
#text!(9, 1, text = L"\times 10 ^{-3}", align = [:right, :bottom], color = :black, fontsize = 30)
lines!(ax, LinRange(1,length(SST),length(SST)), SST.*1e-2, color = :indigo, linewidth = 3)
save("../results/financial/Diks18.png", fig)

