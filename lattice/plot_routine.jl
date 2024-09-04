using CairoMakie
using CSV, DataFrames
using Statistics

CairoMakie.activate!(; px_per_unit = 3)

df = DataFrame(CSV.File("./1/aggregate.csv"; delim=','))
aggregate1 = df[!,1]

start = 1
finish = length(aggregate1) 
μ = LinRange(0.5, 1.2, length(aggregate1))

fig = Figure(; size = (1200, 500), backgroundcolor = :transparent)
ax = Axis(fig[1,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((0.5,1.2), (0.425,0.95)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = true,
                xlabelsize = 36,
                xlabelcolor = :black,
                xlabelpadding = -40.0,
                xticks = [0.5, 1.2],
                xticksvisible = true,
                xticksize = 10,
                xticklabelsvisible = true,
                xticklabelsize = 28,
                xticklabelalign = (:right,:top),
                xtickformat = "{:.1f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\tilde{u}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = -40.0,
                yticks = [0.4,0.425,0.95],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :left,
)
lines!(ax, μ[start:finish], aggregate1[start:finish], linewidth = 2.0, color = (:black, 0.75))
save("./slide10.1.png", fig)

# Detailed region
start = 5775
finish = 5825 
μ = LinRange(0.5, 1.2, length(aggregate1))

fig = Figure(; size = (1200, 500), backgroundcolor = :transparent)
ax = Axis(fig[1,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((μ[start],μ[finish]), (8.1625,8.215)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = true,
                xlabelsize = 36,
                xlabelcolor = :black,
                xlabelpadding = -40.0,
                xticks = [μ[start],μ[finish]],
                xticksvisible = true,
                xticksize = 10,
                xticklabelsvisible = true,
                xticklabelsize = 28,
                xticklabelalign = (:left,:top),
                xtickformat = "{:.2f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\tilde{u}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = -40.0,
                yticks = [8.1625,8.215],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.2f}",
                yscale = identity,
                yaxisposition = :right,
)
lines!(ax, μ[start:finish], aggregate1[start:finish].*10, linewidth = 2.0, color = (:black, 0.75))
save("./slide10.1.bis.png", fig)

# Plot 2
df = DataFrame(CSV.File("./2/aggregate.csv"; delim=','))
aggregate2 = df[!,1]

start = 1
finish = length(aggregate2) 
μ = LinRange(0.5, 1.2, length(aggregate2))

fig = Figure(; size = (1200, 500), backgroundcolor = :transparent)
ax = Axis(fig[1,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((0.5,1.2), (0.425,0.95)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = true,
                xlabelsize = 36,
                xlabelcolor = :black,
                xlabelpadding = -40.0,
                xticks = [0.5, 1.2],
                xticksvisible = true,
                xticksize = 10,
                xticklabelsvisible = true,
                xticklabelsize = 28,
                xticklabelalign = (:right,:top),
                xtickformat = "{:.1f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\tilde{u}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = -40.0,
                yticks = [0.4,0.425,0.95],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :left,
)
lines!(ax, μ[start:finish], aggregate1[start:finish], linewidth = 2.0, color = (:black, 0.75))
lines!(ax, μ[start:finish], aggregate2[start:finish], linewidth = 2.0, color = (:red, 0.5))
save("./slide10.2.png", fig)

# Residuals 
df = DataFrame(CSV.File("./2/imfs4.csv"; delim=','))
residuals2 = df[!,1]
start = 1 
finish = length(aggregate2) 
μ = LinRange(0.5, 1.2, length(aggregate2))

fig = Figure(; size = (1200, 500), backgroundcolor = :transparent)
ax = Axis(fig[1,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((μ[start],μ[finish]), (-1.5,1.5)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = true,
                xlabelsize = 36,
                xlabelcolor = :black,
                xlabelpadding = -40.0,
                xticks = [μ[start],μ[finish]],
                xticksvisible = true,
                xticksize = 10,
                xticklabelsvisible = true,
                xticklabelsize = 28,
                xticklabelalign = (:left,:top),
                xtickformat = "{:.1f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\tilde{u}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = -40.0,
                yticks = [-1.5,1.5],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :right,
)
lines!(ax, μ[start:finish], residuals2[start:finish].*1e5, linewidth = 2.0, color = (:red, 0.5))
text!((1.1,1.15),text=L"\times 10^{-5}",fontsize=36)
save("./slide10.2.bis.png", fig)

# Plot 3
df = DataFrame(CSV.File("./3/aggregate.csv"; delim=','))
aggregate3 = df[!,1]

start = 1
finish = length(aggregate3)
μ = LinRange(0.5, 1.2, length(aggregate3))

fig = Figure(; size = (1200, 500), backgroundcolor = :transparent)
ax = Axis(fig[1,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((0.5,1.2), (0.425,0.95)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = true,
                xlabelsize = 36,
                xlabelcolor = :black,
                xlabelpadding = -40.0,
                xticks = [0.5, 1.2],
                xticksvisible = true,
                xticksize = 10,
                xticklabelsvisible = true,
                xticklabelsize = 28,
                xticklabelalign = (:right,:top),
                xtickformat = "{:.1f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\tilde{u}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = -40.0,
                yticks = [0.4,0.425,0.95],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :left,
)
lines!(ax, μ[start:finish], aggregate1[start:finish], linewidth = 2.0, color = (:black, 0.75))
lines!(ax, μ[start:finish], aggregate2[start:finish], linewidth = 2.0, color = (:red, 0.5))
lines!(ax, μ[start:finish], aggregate3[start:finish], linewidth = 2.0, color = (:blue, 0.5))
save("./slide10.3.png", fig)

# Residuals 
df = DataFrame(CSV.File("./3/imfs7.csv"; delim=','))
residuals3 = df[!,1]
start = 1 
finish = length(aggregate2) 
μ = LinRange(0.5, 1.2, length(aggregate2))

fig = Figure(; size = (1200, 500), backgroundcolor = :transparent)
ax = Axis(fig[1,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((μ[start],μ[finish]), (-1.5,1.5)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = true,
                xlabelsize = 36,
                xlabelcolor = :black,
                xlabelpadding = -40.0,
                xticks = [μ[start],μ[finish]],
                xticksvisible = true,
                xticksize = 10,
                xticklabelsvisible = true,
                xticklabelsize = 28,
                xticklabelalign = (:left,:top),
                xtickformat = "{:.1f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\tilde{u}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = -40.0,
                yticks = [-1.5,1.5],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :right,
)
lines!(ax, μ[start:finish], residuals2[start:finish].*1e4, linewidth = 2.0, color = (:red, 0.5))
lines!(ax, μ[start:finish], residuals3[start:finish].*1e4, linewidth = 2.0, color = (:blue, 0.25))
text!((1.1,1.15),text=L"\times 10^{-4}",fontsize=36)
save("./slide10.3.bis.png", fig)

# Plot 4
df = DataFrame(CSV.File("./4/aggregate.csv"; delim=','))
aggregate4 = df[!,1]

start = 1
finish = length(aggregate4)
μ = LinRange(0.5, 1.2, length(aggregate4))

fig = Figure(; size = (1200, 500), backgroundcolor = :transparent)
ax = Axis(fig[1,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((0.5,1.2), (0.425,0.95)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = true,
                xlabelsize = 36,
                xlabelcolor = :black,
                xlabelpadding = -40.0,
                xticks = [0.5, 1.2],
                xticksvisible = true,
                xticksize = 10,
                xticklabelsvisible = true,
                xticklabelsize = 28,
                xticklabelalign = (:right,:top),
                xtickformat = "{:.1f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\tilde{u}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = -40.0,
                yticks = [0.4,0.425,0.95],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :left,
)
lines!(ax, μ[start:finish], aggregate1[start:finish], linewidth = 2.0, color = (:black, 0.75))
lines!(ax, μ[start:finish], aggregate2[start:finish], linewidth = 2.0, color = (:red, 0.5))
lines!(ax, μ[start:finish], aggregate3[start:finish], linewidth = 2.0, color = (:blue, 0.5))
lines!(ax, μ[start:finish], aggregate4[start:finish], linewidth = 2.0, color = (:green, 0.5))
save("./slide10.4.png", fig)

# Residuals 
df = DataFrame(CSV.File("./4/imfs6.csv"; delim=','))
residuals4 = df[!,1]
start = 1 
finish = length(aggregate2) 
μ = LinRange(0.5, 1.2, length(aggregate2))

fig = Figure(; size = (1200, 500), backgroundcolor = :transparent)
ax = Axis(fig[1,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((μ[start],μ[finish]), (-3.5,3.5)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = true,
                xlabelsize = 36,
                xlabelcolor = :black,
                xlabelpadding = -40.0,
                xticks = [μ[start],μ[finish]],
                xticksvisible = true,
                xticksize = 10,
                xticklabelsvisible = true,
                xticklabelsize = 28,
                xticklabelalign = (:left,:top),
                xtickformat = "{:.1f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\tilde{u}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = -40.0,
                yticks = [-3.5,3.5],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :right,
)
lines!(ax, μ[start:finish], residuals2[start:finish].*1e4, linewidth = 2.0, color = (:red, 0.5))
lines!(ax, μ[start:finish], residuals3[start:finish].*1e4, linewidth = 2.0, color = (:blue, 0.25))
lines!(ax, μ[start:finish], residuals4[start:finish].*1e4, linewidth = 2.0, color = (:green, 0.15))
text!((1.1,2.7),text=L"\times 10^{-4}",fontsize=36)
save("./slide10.4.bis.png", fig)

# Compute the variance of the 3 residuals
width = 200
variance2 = [var(residuals2[n:n+width]) for n in 1:(length(residuals2)-width)]
variance3 = [var(residuals3[n:n+width]) for n in 1:(length(residuals3)-width)]
variance4 = [var(residuals4[n:n+width]) for n in 1:(length(residuals4)-width)]

fig = Figure(; size = (800, 800), backgroundcolor = :transparent)
ax = Axis(fig[1,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((μ[width+1],μ[finish]), nothing),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = true,
                xlabelsize = 36,
                xlabelcolor = :black,
                xlabelpadding = -40.0,
                xticks = [μ[width+1],μ[finish]],
                xticksvisible = true,
                xticksize = 10,
                xticklabelsvisible = true,
                xticklabelsize = 28,
                xticklabelalign = (:right,:top),
                xtickformat = "{:.1f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\text{variance}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = -40.0,
                yticks = [0.0,0.8],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :left,
)
lines!(ax, μ[width+1:finish], variance2.*1e8, linewidth = 2.0, color = (:red, 0.5))
lines!(ax, μ[width+1:finish], variance3.*1e8, linewidth = 2.0, color = (:blue, 0.25))
lines!(ax, μ[width+1:finish], variance4.*1e8, linewidth = 2.0, color = (:green, 0.15))
text!((1.075,0.8),text=L"\times 10^{-8}",fontsize=36)
save("./slide10.5.png", fig)
