using DifferentialEquations
using LaTeXStrings, CairoMakie
using Statistics
using PyCall

#     Ornstein-Uhlenbeck process
# (linear drift, constant diffusion)

# Define the parameters of the process
σ = 0.1 
ε = 0.02
α = 1.0
trueVar = 0.005
∂N = 0.35
x0 = [0.0, -2.0]
T = 10.00
δt = 1e-4

# Define the deterministic dynamics
function iip_det!(f, x, y, t)
        f[1] = -(α/ε)*x[1] 
        f[2] = 1.0 
        return nothing
end
# Define the stochastic dynamics
function iip_stoc!(f, x, y, t)
        f[1] = + σ/sqrt(ε)
        f[2] = 0.0
        return nothing
end

# Build the (non-parametric) dynamical system
OUP = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T))

# Evolve a sample path of the stochastic process 
Xt = solve(OUP, EM(), dt=δt)

# Plot the sample path 
CairoMakie.activate!(; px_per_unit=3)
fig1 = Figure(; size = (1000, 400), backgroundcolor = :transparent)
ax = Axis(fig1[1,2:4],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((-2,8), (-0.5,0.5)),
    # Title
    #title = L"r = r_c = %$R",
    titlevisible = false,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"\mu",
    xlabelvisible = true,
    xlabelsize = 20,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    xticks = [-2,0,2,4,6,8],
    xticksvisible = true,
    xticksize = 6,
    xticklabelsvisible = true,
    xticklabelsize = 14,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"x",
    ylabelvisible = true,
    ylabelsize = 20,
    ylabelcolor = :black,
    ylabelpadding = 10,
    yticks = [-0.5,-0.35,0,0.35,0.5],
    yticksvisible = true,
    yticksize = 6,
    yticklabelsvisible = false,
    yticklabelsize = 14,
    ytickformat = "{:.2f}",
    yscale = identity,
    yaxisposition = :left,
)
limits!(ax,x0[2],x0[2]+T,-∂N-0.15,∂N+0.15)
y_var = range(x0[2],x0[2]+T,length=1000) 
x_var = trueVar .+ 0.0*y_var 
x_∂N = ∂N .+ 0.0*y_var
lines!(y_var, 0.0.*x_var, color = :grey, linewidth = 10)
lines!(y_var, x_∂N, color = :blue, linewidth = 1, linestyle = :dash)
text!(0.0, +0.4, text=L"\partial N^{+}(\rho;C_{\epsilon})", color = :blue, align = (:center, :center))
lines!(y_var, -x_∂N, color = :blue, linewidth = 1, linestyle = :dash)
text!(0.0, -0.4, text=L"\partial N^{-}(\rho;C_{\epsilon})", color = :blue, align = (:center, :center))
lines!(ax, Xt[2,:], Xt[1,:], color= :black, linewidth = 1)

# Plot distribution of the time-series
ax2 = Axis(fig1[1,1],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = (nothing, (-0.5,+0.5)),
    # Title
    #title = L"r = r_c = %$R",
    titlevisible = false,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"\text{density}",
    xlabelvisible = true,
    xlabelsize = 20,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    xticks = [0,0.06],
    xticksvisible = true,
    xticksize = 6,
    xticklabelsvisible = true,
    xticklabelsize = 14,
    xtickformat = "{:.2f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"x",
    ylabelvisible = false,
    ylabelsize = 20,
    ylabelcolor = :black,
    ylabelpadding = -15.0,
    yticks = [-0.5,-0.35,0,0.35,0.5],
    yticksvisible = true,
    yticksize = 6,
    yticklabelsvisible = true,
    yticklabelsize = 14,
    ytickformat = "{:.2f}",
    yscale = identity,
    yaxisposition = :left,
)
hist!(ax2, Xt[1,:], bins = 50, normalization = :probability, color = :red, strokecolor = :black, strokewidth = 1, direction = :x)

# Plot the time-varying moments of the distribution 
fig2 = Figure(; size = (1000, 400), backgroundcolor = :transparent)
ax3 = Axis(fig2[1,1],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((-2,8), (0,0.01)),
    # Title
    #title = L"r = r_c = %$R",
    titlevisible = false,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"\mu",
    xlabelvisible = true,
    xlabelsize = 20,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    xticks = [-2,0,2,4,6,8],
    xticksvisible = true,
    xticksize = 6,
    xticklabelsvisible = true,
    xticklabelsize = 14,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"\text{variance}",
    ylabelvisible = true,
    ylabelsize = 20,
    ylabelcolor = :black,
    ylabelpadding = 5,
    yticks = [0, 0.005, 0.01],
    yticksvisible = true,
    yticksize = 6,
    yticklabelsvisible = true,
    yticklabelsize = 14,
    ytickformat = "{:.3f}",
    yscale = identity,
    yaxisposition = :left,
)
time = LinRange(0,T,size(Xt,2)) 
# Numerical values computed from time-series
E = mean(Xt[1,:]) .+ 0.0*time
Var = zeros(Float64, size(Xt,2))
for n in 1:size(Xt,2)
        Var[n] = var(Xt[1,1:n]; mean = E[1])
end
# Theoretical values [Kuehn, 2011]
E_exact = x0[1].*exp.((-α.*time)./ε)
Var_exact = (x0[1] - (σ^2)/(2*α)).*exp.((-2*α.*time)./ε) .+ (σ^2)/(2*α)
l1 = lines!(ax3, Xt[2,:], Var, color = :red)
l2 = lines!(ax3, Xt[2,:], Var_exact, color = :blue, linestyle = :dash)

# Export the figure
save("../../results/fast_slow/fig2.3.png", fig1)
save("../../results/fast_slow/fig2.4.png", fig2)
