using DynamicalSystems, DifferentialEquations
using LaTeXStrings, CairoMakie, Makie.Colors
using Statistics, ProgressMeter 

# Noise level and timescale separation
σ = 0.05
ε = 0.01
# Timestep
δt = 1e-1
T = 400.00
# IC 
x0 = [3.06,-1.0]

# Deterministic dynamics 
function iip_det!(f, x, y, t)
        f[1] = -x[2] - 2*x[1] + 3*(x[1])^2 - 0.8*(x[1])^3
        f[2] = ε
        return nothing
end
# Stochastic dynamics
function iip_stoc!(f, x, y, t)
        f[1] = +σ
        f[2] = 0.0
        return nothing
end

# Define and solve the problem 
normal_form = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T))
sol = solve(normal_form, EM(), dt=δt, verbose=false)
xt = sol[1,:]
μt = sol[2,:]

# Define boundaries for plotting purposes
xt_max = xt[argmax(xt)] + 0.5
xt_min = xt[argmin(xt)] - 0.5

# Plot the timeseries
CairoMakie.activate!(; px_per_unit = 3)
fig1 = Figure(; size = (1200, 600), backgroundcolor = :transparent)
ax1 = Axis(fig1[1,1:4],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((0, T*1e-3), (xt_min, xt_max)),
    # Title
    title = "Title",
    titlevisible = false,
    titlesize = 25,
    titlealign = :center,
    titlegap = -38.0,
    # x-axis
    xlabel = L"\text{time}",
    xlabelvisible = true,
    xlabelsize = 30,
    xlabelcolor = :black,
    xlabelpadding = 10.0,
    #xticks = [-3,-2,-1,0],
    xticksvisible = false,
    xticksize = 20,
    xticklabelsvisible = false,
    xticklabelsize = 20,
    #xtickformat = values -> [L"\sqrt{%$(round(value; digits = 1))}" for value in values],
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"\text{observable}",
    ylabelvisible = true,
    ylabelsize = 30,
    ylabelcolor = :black,
    ylabelpadding = 10.0,
    #yticks = [-0.5,0.5],
    yticksvisible = false,
    yticksize = 20,
    yticklabelsvisible = false,
    yticklabelsize = 20,
    ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
)
lines!(ax1, LinRange(0.0, T, length(xt)).*1e-3, xt, linewidth = 1.0, color = (:red, 0.5))
save("../../results/precursors/slide2.1.png", fig1)

# Plot the state-space trajectory
CairoMakie.activate!(; px_per_unit = 3)
fig2 = Figure(; size = (1200, 600), backgroundcolor = :transparent)
ax2 = Axis(fig2[1,1:4],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((μt[1], μt[end]), (xt_min, xt_max)),
    # Title
    title = "Title",
    titlevisible = false,
    titlesize = 25,
    titlealign = :center,
    titlegap = -38.0,
    # x-axis
    xlabel = L"\mu",
    xlabelvisible = true,
    xlabelsize = 30,
    xlabelcolor = :black,
    xlabelpadding = 10.0,
    #xticks = [-3,-2,-1,0],
    xticksvisible = false,
    xticksize = 20,
    xticklabelsvisible = false,
    xticklabelsize = 20,
    #xtickformat = values -> [L"\sqrt{%$(round(value; digits = 1))}" for value in values],
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"x",
    ylabelvisible = true,
    ylabelsize = 30,
    ylabelcolor = :black,
    ylabelpadding = 10.0,
    #yticks = [-0.5,0.5],
    yticksvisible = false,
    yticksize = 20,
    yticklabelsvisible = false,
    yticklabelsize = 20,
    ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
)
# Plot the bifurcation diagram
a = μt[1]
b = μt[end]
μ=LinRange(a,b,1000)
x1 = Float64[]
μ1 = Float64[]
x2 = Float64[]
μ2 = Float64[]
x3 = Float64[]
μ3 = Float64[]
equilibria = ComplexF64[3]
using Polynomials
for n in 1:length(μ) 
        global equilibria = roots(Polynomial([-μ[n],-2,+3,-0.8]))
        if imag(equilibria[1])≈0.0
                push!(x1, real(equilibria[1]))
                push!(μ1, μ[n])
        end
        if imag(equilibria[2])≈0.0
                push!(x2, real(equilibria[2]))
                push!(μ2, μ[n])
        end
        if imag(equilibria[3])≈0.0
                push!(x3, real(equilibria[3]))
                push!(μ3, μ[n])
        end
end
lines!(ax2, μ1, x1, color = :black, linewidth = 1.5)
lines!(ax2, μ2, x2, color = :black, linestyle = :dash, linewidth = 1.5)
lines!(ax2, μ3, x3, color = :black, linewidth = 1.5)
lines!(ax2, μt, xt, linewidth = 1.0, color = (:red, 0.5))
save("../../results/precursors/slide2.2.1.png", fig2)

# Plot the variance 
width = 400
variance = [var(xt[n:n+width]) for n in 1:(length(xt)-width)]
c = -0.005
d = 0.05

CairoMakie.activate!(; px_per_unit = 3)
fig3 = Figure(; size = (1200, 600), backgroundcolor = :transparent)
ax3 = Axis(fig3[1,1:4],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((a, b), (c, d)),
    # Title
    title = "Title",
    titlevisible = false,
    titlesize = 25,
    titlealign = :center,
    titlegap = -38.0,
    # x-axis
    xlabel = L"\text{time}",
    xlabelvisible = true,
    xlabelsize = 30,
    xlabelcolor = :black,
    xlabelpadding = 10.0,
    #xticks = [-3,-2,-1,0],
    xticksvisible = false,
    xticksize = 20,
    xticklabelsvisible = false,
    xticklabelsize = 20,
    #xtickformat = values -> [L"\sqrt{%$(round(value; digits = 1))}" for value in values],
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"\text{variance}",
    ylabelvisible = true,
    ylabelsize = 30,
    ylabelcolor = :black,
    ylabelpadding = 10.0,
    #yticks = [-0.5,0.5],
    yticksvisible = false,
    yticksize = 20,
    yticklabelsvisible = false,
    yticklabelsize = 20,
    #ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
)
lines!(ax3, 1.62.*ones(2), [c, d], color = :black, linewidth = 2, linestyle = :dash)
lines!(ax3, LinRange(μt[width],b,length(variance)), variance, color = :blue, linewidth = 1.25)
save("../../results/precursors/slide2.3.1.png", fig3)

# Plot the initial position of the sliding window for both the timeseries and variance
using Makie.GeometryBasics
fig4 = Figure(; size = (1200, 600), backgroundcolor = :transparent)
ax4 = Axis(fig4[1,1:4],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((a, b), (c, d)),
    # Title
    title = "Title",
    titlevisible = false,
    titlesize = 25,
    titlealign = :center,
    titlegap = -38.0,
    # x-axis
    xlabel = L"\text{time}",
    xlabelvisible = true,
    xlabelsize = 30,
    xlabelcolor = :black,
    xlabelpadding = 10.0,
    #xticks = [-3,-2,-1,0],
    xticksvisible = false,
    xticksize = 20,
    xticklabelsvisible = false,
    xticklabelsize = 20,
    #xtickformat = values -> [L"\sqrt{%$(round(value; digits = 1))}" for value in values],
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"\text{variance}",
    ylabelvisible = true,
    ylabelsize = 30,
    ylabelcolor = :black,
    ylabelpadding = 10.0,
    #yticks = [-0.5,0.5],
    yticksvisible = false,
    yticksize = 20,
    yticklabelsvisible = false,
    yticklabelsize = 20,
    #ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
)
lines!(ax4, 1.62.*ones(2), [c, d], color = :black, linewidth = 2, linestyle = :dash)
poly!(ax4, Point2f[(μt[1], c), (μt[width], c), (μt[width], d), (μt[1], d)], color = (:grey, 0.25), strokecolor = :grey, strokewidth = 0.05)
scatter!(ax4, μt[width], variance[1], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 15)
save("../../results/precursors/slide2.3.2.png", fig4)

fig5 = Figure(; size = (1200, 600), backgroundcolor = :transparent)
ax5 = Axis(fig5[1,1:4],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((μt[1], μt[end]), (xt_min, xt_max)),
    # Title
    title = "Title",
    titlevisible = false,
    titlesize = 25,
    titlealign = :center,
    titlegap = -38.0,
    # x-axis
    xlabel = L"\text{time}",
    xlabelvisible = true,
    xlabelsize = 30,
    xlabelcolor = :black,
    xlabelpadding = 10.0,
    #xticks = [-3,-2,-1,0],
    xticksvisible = false,
    xticksize = 20,
    xticklabelsvisible = false,
    xticklabelsize = 20,
    #xtickformat = values -> [L"\sqrt{%$(round(value; digits = 1))}" for value in values],
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"\text{observable}",
    ylabelvisible = true,
    ylabelsize = 30,
    ylabelcolor = :black,
    ylabelpadding = 10.0,
    #yticks = [-0.5,0.5],
    yticksvisible = false,
    yticksize = 20,
    yticklabelsvisible = false,
    yticklabelsize = 20,
    ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
)
# Plot the bifurcation diagram
lines!(ax5, μ1, x1, color = :black, linewidth = 1.5)
lines!(ax5, μ2, x2, color = :black, linestyle = :dash, linewidth = 1.5)
lines!(ax5, μ3, x3, color = :black, linewidth = 1.5)
lines!(ax5, μt, xt, linewidth = 1.0, color = (:red, 0.5))
poly!(ax5, Point2f[(μt[1], xt_min), (μt[width], xt_min), (μt[width], xt_max), (μt[1], xt_max)], color = (:grey, 0.25), strokecolor = :grey, strokewidth = 0.05)
save("../../results/precursors/slide2.2.2.png", fig5)

# Make figures for animation of the sliding window
CairoMakie.activate!(; px_per_unit = 1)
for n in 1:length(variance)
        fig6 = Figure(; size = (1200, 600), backgroundcolor = "#fdfff2ff")
        ax6 = Axis(fig6[1,1:4],
                # Background
                backgroundcolor = "#fdfff2ff",
                xgridvisible = false,
                ygridvisible = false,
                limits = ((μt[1], μt[end]), (xt_min, xt_max)),
                # Title
                title = "Title",
                titlevisible = false,
                titlesize = 25,
                titlealign = :center,
                titlegap = -38.0,
                # x-axis
                xlabel = L"\text{time}",
                xlabelvisible = true,
                xlabelsize = 30,
                xlabelcolor = :black,
                xlabelpadding = 10.0,
                #xticks = [-3,-2,-1,0],
                xticksvisible = false,
                xticksize = 20,
                xticklabelsvisible = false,
                xticklabelsize = 20,
                #xtickformat = values -> [L"\sqrt{%$(round(value; digits = 1))}" for value in values],
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\text{observable}",
                ylabelvisible = true,
                ylabelsize = 30,
                ylabelcolor = :black,
                ylabelpadding = 10.0,
                #yticks = [-0.5,0.5],
                yticksvisible = false,
                yticksize = 20,
                yticklabelsvisible = false,
                yticklabelsize = 20,
                ytickformat = "{:.0f}",
                yscale = identity,
                yaxisposition = :left,
        )
        # Plot the bifurcation diagram
        lines!(ax6, μ1, x1, color = :black, linewidth = 1.5)
        lines!(ax6, μ2, x2, color = :black, linestyle = :dash, linewidth = 1.5)
        lines!(ax6, μ3, x3, color = :black, linewidth = 1.5)
        lines!(ax6, μt, xt, linewidth = 1.0, color = (:red, 0.5))
        poly!(ax6, Point2f[(μt[n], xt_min), (μt[width+n], xt_min), (μt[width+n], xt_max), (μt[n], xt_max)], color = (:grey, 0.25), strokecolor = :grey, strokewidth = 0.05)
        save("../../results/precursors/animation/timeseries/$n.png", fig6)

        # Variance
        fig7 = Figure(; size = (1200, 600), backgroundcolor = "#fdfff2ff")
        ax7 = Axis(fig7[1,1:4],
                # Background
                backgroundcolor = "#fdfff2ff",
                xgridvisible = false,
                ygridvisible = false,
                limits = ((a, b), (c, d)),
                # Title
                title = "Title",
                titlevisible = false,
                titlesize = 25,
                titlealign = :center,
                titlegap = -38.0,
                # x-axis
                xlabel = L"\text{time}",
                xlabelvisible = true,
                xlabelsize = 30,
                xlabelcolor = :black,
                xlabelpadding = 10.0,
                #xticks = [-3,-2,-1,0],
                xticksvisible = false,
                xticksize = 20,
                xticklabelsvisible = false,
                xticklabelsize = 20,
                #xtickformat = values -> [L"\sqrt{%$(round(value; digits = 1))}" for value in values],
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\text{variance}",
                ylabelvisible = true,
                ylabelsize = 30,
                ylabelcolor = :black,
                ylabelpadding = 10.0,
                #yticks = [-0.5,0.5],
                yticksvisible = false,
                yticksize = 20,
                yticklabelsvisible = false,
                yticklabelsize = 20,
                ytickformat = "{:.0f}",
                yscale = identity,
                yaxisposition = :left,
        )
        lines!(ax7, 1.62.*ones(2), [c, d], color = :black, linewidth = 2, linestyle = :dash)
        lines!(ax7, LinRange(μt[width],μt[width+n],n+1), variance[1:n+1], color = :blue, linewidth = 1.25)
        poly!(ax7, Point2f[(μt[n], c), (μt[width+n], c), (μt[width+n], d), (μt[n], d)], color = (:grey, 0.25), strokecolor = :grey, strokewidth = 0.05)
        scatter!(ax7, μt[width+n], variance[n], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 15)
        save("../../results/precursors/animation/variance/$n.png", fig7)
end
