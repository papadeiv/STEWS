using CairoMakie
using CSV, DataFrames
using Statistics

CairoMakie.activate!(; px_per_unit = 3)

df = DataFrame(CSV.File("./nodetrend.csv"; delim=','))
variance = df[!,1]
fig = Figure(; size = (1200, 600), backgroundcolor = :transparent)
ax = Axis(fig[1,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                #limits = ((0.5,1.2), (1.28,3.25)),
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
                xlabelpadding = -20.0,
                #xticks = [0.5, 0.68, 1.2],
                xticksvisible = true,
                xticksize = 10,
                xticklabelsvisible = true,
                xticklabelsize = 28,
                xticklabelalign = (:right,:top),
                #xtickformat = "{:.2f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\text{variance}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = -30.0,
                #yticks = [1.3, 3.1],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                #ytickformat = "{:.2f}",
                yscale = identity,
                yaxisposition = :left,
)
lines!(ax, LinRange(0.57,0.68,length(variance)), variance, linewidth = 2.0)
#lines!(ax, 0.68.*ones(2), [1.28,3.25], color = :black, linestyle = :dash)
#text!((0.525,3.00),text=L"\times 10^{-5}",fontsize=46)
save("./slide5.1.png", fig)

df = DataFrame(CSV.File("./imfs1.csv"; delim=','))
imfs1 = df[!,1]
fig = Figure(; size = (1200, 600), backgroundcolor = :transparent)#"#fdfff2ff")
ax = Axis(fig[1,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                #limits = ((0.5,1.2), (0.5,1.0)),
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
                xlabelpadding = -20.0,
                #xticks = [0.5, 0.68, 1.2],
                xticksvisible = true,
                xticksize = 15,
                xticklabelsvisible = true,
                xticklabelsize = 28,
                xticklabelalign = (:right,:top),
                #xtickformat = "{:.2f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\text{trend}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = -30.0,
                #yticks = [0.5,1.0],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                #ytickformat = "{:.1f}",
                yscale = identity,
                yaxisposition = :left,
)
lines!(ax, LinRange(0.5,1.2,length(imfs1)), imfs1, linewidth = 2.0)
#lines!(ax, 0.68.*ones(2), [0.5,1], color = :black, linestyle = :dash)
save("./slide5.2.png", fig)

#=
df = DataFrame(CSV.File("./imfs2.csv"; delim=','))
imfs2 = df[!,1]
fig = Figure(; size = (1200, 600))#, backgroundcolor = :transparent)#"#fdfff2ff")
ax = Axis(fig[1,1],
           #title = L"\mu = %$parameter",
           titlesize = 25,
           #backgroundcolor = :transparent,#"#fdfff2ff"
           xgridvisible = false,
           ygridvisible = false,
           #xticksvisible = false,
           #yticksvisible = false,
           #xticklabelsvisible = false,
           #yticklabelsvisible = false,
           limits = (0.5, 1.2, nothing, nothing)
)
lines!(ax, LinRange(0.5,1.2,length(imfs2)), imfs2, linewidth = 2.0)
save("../results/lattice/lung_ventilation/imfs2.png", fig)

width = floor(Int, 0.1*length(imfs2))
variance_imfs2 = [var(imfs2[t:t+width]) for t in 1:(length(imfs2)-width)]
fig = Figure(; size = (1200, 600))#, backgroundcolor = :transparent)#"#fdfff2ff")
ax = Axis(fig[1,1],
           #title = L"\mu = %$parameter",
           titlesize = 25,
           #backgroundcolor = :transparent,#"#fdfff2ff"
           xgridvisible = false,
           ygridvisible = false,
           limits = (0.5, 1.2, nothing, nothing)
)
lines!(ax, LinRange(0.57,1.2,length(variance_imfs2)), variance_imfs2, linewidth = 2.0)
lines!(ax, 0.68.*ones(2), [variance_imfs2[argmin(variance_imfs2)], variance_imfs2[argmax(variance_imfs2)]], color = :black, linestyle = :dash)
Label(fig[1, 1, Top()], halign = :left, L"\times 10^{-7}")
save("../results/lattice/lung_ventilation/var_imfs2.png", fig)

df = DataFrame(CSV.File("./imfs3.csv"; delim=','))
imfs3 = df[!,1]
fig = Figure(; size = (1200, 600))#, backgroundcolor = :transparent)#"#fdfff2ff")
ax = Axis(fig[1,1],
           #title = L"\mu = %$parameter",
           titlesize = 25,
           #backgroundcolor = :transparent,#"#fdfff2ff"
           xgridvisible = false,
           ygridvisible = false,
           limits = (0.5, 1.2, nothing, nothing)
)
lines!(ax, LinRange(0.5,1.2,length(imfs3)), imfs3, linewidth = 2.0)
save("../results/lattice/lung_ventilation/imfs3.png", fig)

width = floor(Int, 0.1*length(imfs3))
variance_imfs3 = [var(imfs3[t:t+width]) for t in 1:(length(imfs3)-width)]
fig = Figure(; size = (1200, 600))#, backgroundcolor = :transparent)#"#fdfff2ff")
ax = Axis(fig[1,1],
           #title = L"\mu = %$parameter",
           titlesize = 25,
           #backgroundcolor = :transparent,#"#fdfff2ff"
           xgridvisible = false,
           ygridvisible = false,
           limits = (0.5, 1.2, nothing, nothing)
)
lines!(ax, LinRange(0.57,1.2,length(variance_imfs3)), variance_imfs3, linewidth = 2.0)
lines!(ax, 0.68.*ones(2), [variance_imfs3[argmin(variance_imfs3)], variance_imfs3[argmax(variance_imfs3)]], color = :black, linestyle = :dash)
Label(fig[1, 1, Top()], halign = :left, L"\times 10^{-7}")
save("../results/lattice/lung_ventilation/var_imfs3.png", fig)

df = DataFrame(CSV.File("./imfs4.csv"; delim=','))
imfs4 = df[!,1]
fig = Figure(; size = (1200, 600))#, backgroundcolor = :transparent)#"#fdfff2ff")
ax = Axis(fig[1,1],
           #title = L"\mu = %$parameter",
           titlesize = 25,
           #backgroundcolor = :transparent,#"#fdfff2ff"
           xgridvisible = false,
           ygridvisible = false,
           limits = (0.5, 1.2, nothing, nothing)
)
lines!(ax, LinRange(0.5,1.2,length(imfs4)), imfs4, linewidth = 2.0)
save("../results/lattice/lung_ventilation/imfs4.png", fig)

width = floor(Int, 0.1*length(imfs4))
variance_imfs4 = [var(imfs4[t:t+width]) for t in 1:(length(imfs4)-width)]
fig = Figure(; size = (1200, 600))#, backgroundcolor = :transparent)#"#fdfff2ff")
ax = Axis(fig[1,1],
           #title = L"\mu = %$parameter",
           titlesize = 25,
           #backgroundcolor = :transparent,#"#fdfff2ff"
           xgridvisible = false,
           ygridvisible = false,
           limits = (0.5, 1.2, nothing, nothing)
)
lines!(ax, LinRange(0.57,1.2,length(variance_imfs4)), variance_imfs4, linewidth = 2.0)
lines!(ax, 0.68.*ones(2), [variance_imfs4[argmin(variance_imfs4)], variance_imfs4[argmax(variance_imfs4)]], color = :black, linestyle = :dash)
Label(fig[1, 1, Top()], halign = :left, L"\times 10^{-6}")
save("../results/lattice/lung_ventilation/var_imfs4.png", fig)

df = DataFrame(CSV.File("./imfs5.csv"; delim=','))
imfs5 = df[!,1]
fig = Figure(; size = (1200, 600))#, backgroundcolor = :transparent)#"#fdfff2ff")
ax = Axis(fig[1,1],
           #title = L"\mu = %$parameter",
           titlesize = 25,
           #backgroundcolor = :transparent,#"#fdfff2ff"
           xgridvisible = false,
           ygridvisible = false,
           limits = (0.5, 1.2, nothing, nothing)
)
lines!(ax, LinRange(0.5,1.2,length(imfs5)), imfs5, linewidth = 2.0)
save("../results/lattice/lung_ventilation/imfs5.png", fig)

width = floor(Int, 0.1*length(imfs5))
variance_imfs5 = [var(imfs5[t:t+width]) for t in 1:(length(imfs5)-width)]
fig = Figure(; size = (1200, 600))#, backgroundcolor = :transparent)#"#fdfff2ff")
ax = Axis(fig[1,1],
           #title = L"\mu = %$parameter",
           titlesize = 25,
           #backgroundcolor = :transparent,#"#fdfff2ff"
           xgridvisible = false,
           ygridvisible = false,
           limits = (0.5, 1.2, nothing, nothing)
)
lines!(ax, LinRange(0.57,1.2,length(variance_imfs5)), variance_imfs5, linewidth = 2.0)
lines!(ax, 0.68.*ones(2), [variance_imfs5[argmin(variance_imfs5)], variance_imfs5[argmax(variance_imfs5)]], color = :black, linestyle = :dash)
Label(fig[1, 1, Top()], halign = :left, L"\times 10^{-7}")
save("../results/lattice/lung_ventilation/var_imfs5.png", fig)

#=
df = DataFrame(CSV.File("./emd.csv"; delim=','))
detrended_variance = df[!,1]
fig = Figure(; size = (1200, 600))#, backgroundcolor = :transparent)#"#fdfff2ff")
ax = Axis(fig[1,1],
           #title = L"\mu = %$parameter",
           titlesize = 25,
           #backgroundcolor = :transparent,#"#fdfff2ff"
           xgridvisible = false,
           ygridvisible = false,
)
lines!(ax, LinRange(0.5,1.2,length(detrended_variance)), detrended_variance, linewidth = 2.0)
save("../results/lattice/lung_ventilation/fig3.2.2.png", fig)
=#
=#

df = DataFrame(CSV.File("./aggregate.csv"; delim=','))
aggregate = df[!,1]
detrended_signal = aggregate - imfs1 
fig = Figure(; size = (1200, 600), backgroundcolor = :transparent)#"#fdfff2ff")
ax = Axis(fig[1,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                #limits = ((0.5,1.2), (-3.5,3.5)),
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
                xlabelpadding = -20.0,
                #xticks = [0.5, 0.68, 1.2],
                xticksvisible = true,
                xticksize = 15,
                xticklabelsvisible = true,
                xticklabelsize = 28,
                xticklabelalign = (:right,:top),
                #xtickformat = "{:.2f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\text{residuals}",
                ylabelvisible = true,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = -30.0,
                #yticks = [-3,3],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 28,
                #ytickformat = "{:.2f}",
                yscale = identity,
                yaxisposition = :left,
)
lines!(ax, LinRange(0.5,1.2,length(detrended_signal)), detrended_signal, linewidth = 2.0)
#lines!(ax, 0.68.*ones(2), [-3,3]*1e3, color = :black, linestyle = :dash)
#text!((0.525,2.75),text=L"\times 10^{-3}",fontsize=46)
save("./slide5.3.png", fig)

width = floor(Int, 0.1*length(detrended_signal))
detrended_variance = [var(detrended_signal[t:t+width]) for t in 1000:(length(detrended_signal)-width)]
fig = Figure(; size = (1200, 400), backgroundcolor = :transparent)
ax = Axis(fig[1,1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                #limits = ((0.5,1.2), (0.03,3.25)), # ((0.5,1.2), (3.6,4.4)),
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
                xlabelpadding = -20.0,
                #xticks = [0.5, 0.68, 1.2],
                xticksvisible = false,
                xticksize = 15,
                xticklabelsvisible = false,
                xticklabelsize = 28,
                xticklabelalign = (:right,:top),
                #xtickformat = "{:.2f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\text{variance}",
                ylabelvisible = false,
                ylabelsize = 36,
                ylabelcolor = :black,
                ylabelpadding = -30.0,
                #yticks = [1.3, 3.1], # [3.7,4.3],
                yticksvisible = false,
                yticksize = 10,
                yticklabelsvisible = false,
                yticklabelsize = 28,
                #ytickformat = "{:.5f}",
                yscale = identity,
                yaxisposition = :left,

)
lines!(ax, LinRange(0.57,0.68,length(detrended_variance)), detrended_variance, linewidth = 4.0, color = (:black, 0.5))
hidespines!(ax)
#lines!(ax, 0.68.*ones(2), [1.28,3.25], color = :black, linestyle = :dash)
#text!((0.525,3.00),text=L"\times 10^{-5}",fontsize=46)
save("./slide5.4.png", fig)
