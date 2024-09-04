using CairoMakie
using CSV, DataFrames
using Statistics

CairoMakie.activate!(; px_per_unit = 3)
fig = Figure(; size = (1800, 1200), backgroundcolor = :transparent)

threshold = 3.04 

df = DataFrame(CSV.File("./imfs7.csv"; delim=','))
imfs1 = df[!,1]
ax1 = Axis(fig[1,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((3,3.1), (-2,2)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = false,
                xlabelsize = 36,
                xlabelcolor = :black,
                xlabelpadding = 0.0,
                xticks = [3, 3.02, 3.04, 3.06, 3.08, 3.1],
                xticksvisible = true,
                xticksize = 15,
                xticklabelsvisible = false,
                xticklabelsize = 28,
                xticklabelalign = (:center,:top),
                xtickformat = "{:.2f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\text{IMF1}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = 0.0,
                yticks = [-2,0,2],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :left,
)
lines!(ax1, LinRange(3,3.1,length(imfs1)), imfs1.*1e6, linewidth = 2.0)
lines!(ax1, threshold.*ones(2), [-2,2], color = :black, linestyle = :dash)
text!(3.015, 1.75, text = L"\times 10^{-6}", align = [:right, :top], color = :black, fontsize = 32)

df = DataFrame(CSV.File("./imfs6.csv"; delim=','))
imfs2 = df[!,1]
ax2 = Axis(fig[2,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((3,3.1), (-3,3)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = false,
                xlabelsize = 36,
                xlabelcolor = :black,
                xlabelpadding = 0.0,
                xticks = [3, 3.02, 3.04, 3.06, 3.08, 3.1],
                xticksvisible = true,
                xticksize = 15,
                xticklabelsvisible = false,
                xticklabelsize = 28,
                xticklabelalign = (:center,:top),
                xtickformat = "{:.2f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\text{IMF2}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = 0.0,
                yticks = [-3,0,3],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :left,
)
lines!(ax2, LinRange(3,3.1,length(imfs2)), imfs2.*1e6, linewidth = 2.0)
lines!(ax2, threshold.*ones(2), [-3,3], color = :black, linestyle = :dash)
text!(3.015, 2.75, text = L"\times 10^{-6}", align = [:right, :top], color = :black, fontsize = 32)

df = DataFrame(CSV.File("./imfs5.csv"; delim=','))
imfs3 = df[!,1]
ax3 = Axis(fig[3,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((3,3.1), (-2,2)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = false,
                xlabelsize = 36,
                xlabelcolor = :black,
                xlabelpadding = 0.0,
                xticks = [3, 3.02, 3.04, 3.06, 3.08, 3.1],
                xticksvisible = true,
                xticksize = 15,
                xticklabelsvisible = false,
                xticklabelsize = 28,
                xticklabelalign = (:center,:top),
                xtickformat = "{:.2f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\text{IMF3}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = 0.0,
                yticks = [-2,0,2],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :left,
)
lines!(ax3, LinRange(3,3.1,length(imfs3)), imfs3.*1e6, linewidth = 2.0)
lines!(ax3, threshold.*ones(2), [-2,2], color = :black, linestyle = :dash)
text!(3.015, 1.75, text = L"\times 10^{-6}", align = [:right, :top], color = :black, fontsize = 32)

df = DataFrame(CSV.File("./imfs4.csv"; delim=','))
imfs4 = df[!,1]
ax4 = Axis(fig[4,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((3,3.1), (-2,2)),
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
                xlabelpadding = 0,
                xticks = [3, 3.02, 3.04, 3.06, 3.08, 3.1],
                xticksvisible = true,
                xticksize = 15,
                xticklabelsvisible = true,
                xticklabelsize = 28,
                xticklabelalign = (:right,:top),
                xtickformat = "{:.2f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\text{IMF4}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = 0.0,
                yticks = [-2,0,2],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :left,
)
lines!(ax4, LinRange(3,3.1,length(imfs4)), imfs4.*1e6, linewidth = 2.0)
lines!(ax4, threshold.*ones(2), [-2,2], color = :black, linestyle = :dash)
text!(3.015, 1.75, text = L"\times 10^{-6}", align = [:right, :top], color = :black, fontsize = 32)

df = DataFrame(CSV.File("./imfs3.csv"; delim=','))
imfs5 = df[!,1]
ax5 = Axis(fig[2,2],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((3,3.1), (-1.25,1.25)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = false,
                xlabelsize = 36,
                xlabelcolor = :black,
                xlabelpadding = 0.0,
                xticks = [3, 3.02, 3.04, 3.06, 3.08, 3.1],
                xticksvisible = true,
                xticksize = 15,
                xticklabelsvisible = false,
                xticklabelsize = 28,
                xticklabelalign = (:center,:top),
                xtickformat = "{:.2f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\text{IMF5}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = 0.0,
                yticks = [-1.25,0,1.25],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :right,
)
lines!(ax5, LinRange(3,3.1,length(imfs5)), imfs5.*1e5, linewidth = 2.0)
lines!(ax5, threshold.*ones(2), [-1.25,1.25], color = :black, linestyle = :dash)
text!(3.0975, 1.1, text = L"\times 10^{-5}", align = [:right, :top], color = :black, fontsize = 32)

df = DataFrame(CSV.File("./imfs2.csv"; delim=','))
imfs6 = df[!,1]
ax6 = Axis(fig[3,2],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((3,3.1), (-4.5,4.5)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = false,
                xlabelsize = 36,
                xlabelcolor = :black,
                xlabelpadding = 0.0,
                xticks = [3, 3.02, 3.04, 3.06, 3.08, 3.1],
                xticksvisible = true,
                xticksize = 15,
                xticklabelsvisible = false,
                xticklabelsize = 28,
                xticklabelalign = (:center,:top),
                xtickformat = "{:.2f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\text{IMF6}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = 0.0,
                yticks = [-4.5,0,4.5],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :right,
)
lines!(ax6, LinRange(3,3.1,length(imfs6)), imfs6.*1e5, linewidth = 2.0)
lines!(ax6, threshold.*ones(2), [-4.5,4.5], color = :black, linestyle = :dash)
text!(3.0975, 3.8, text = L"\times 10^{-5}", align = [:right, :top], color = :black, fontsize = 32)

df = DataFrame(CSV.File("./imfs1.csv"; delim=','))
imfs7 = df[!,1]
ax7 = Axis(fig[4,2],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((3,3.1), (0.75,1.9)),
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
                xlabelpadding = 0,
                xticks = [3, 3.02, 3.04, 3.06, 3.08, 3.1],
                xticksvisible = true,
                xticksize = 15,
                xticklabelsvisible = true,
                xticklabelsize = 28,
                xticklabelalign = (:left,:top),
                xtickformat = "{:.2f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\text{IMF7}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = 0.0,
                yticks = [0.75,1.9],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :right,
)
lines!(ax7, LinRange(3,3.1,length(imfs7)), imfs7.*1e4, linewidth = 2.0)
lines!(ax7, threshold.*ones(2), [0.75,1.9], color = :black, linestyle = :dash)
text!(3.0975, 1.8, text = L"\times 10^{-4}", align = [:right, :top], color = :black, fontsize = 32)

save("../../results/lattice/lung_ventilation/fig2.9.png", fig)
