using DynamicalSystems, DifferentialEquations
using LaTeXStrings, CairoMakie, Makie.Colors
using Statistics, ProgressMeter 

# Noise level and timescale separation
σ = 1.25
ε = 0.01
# Initial state
x0 = [3.33,-3.0]
# Temporal evolution
T = 1000.00
δt = 1e-3

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
# Definition of the ensemble problem
normal_form = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T))
Ne = 10

# Arrays to store the ensemble sample paths
xt = fill(Float64[], 0) 
μt = fill(Float64[], 0)

# Solution of the ensemble problem 
println("Simulating the ensemble sample paths")
@showprogress for j in 1:Ne
        local sol = solve(normal_form, EM(), dt=δt, verbose=false)
        push!(xt, sol[1,:])
        push!(μt, sol[2,:])
end

# Plot the solutions
CairoMakie.activate!(; px_per_unit = 3)
fig = Figure(; size = (1200, 800), backgroundcolor = :transparent)
ax = Axis(fig[1:4,1:3],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((-3,4), nothing),
    # Title
    title = "Title",
    titlevisible = false,
    titlesize = 25,
    titlealign = :center,
    titlegap = -38.0,
    # x-axis
    xlabel = L"\mu",
    xlabelvisible = true,
    xlabelsize = 38,
    xlabelcolor = :black,
    xlabelpadding = -42.0,
    xticks = [-3,-2,-1,0,1,2,3,4],
    xticksvisible = true,
    xticksize = 20,
    xticklabelsvisible = true,
    xticklabelsize = 34,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"x",
    ylabelvisible = true,
    ylabelsize = 38,
    ylabelcolor = :black,
    ylabelpadding = -38.0,
    yticks = [-2,-1,0,2,3,4],
    yticksvisible = true,
    yticksize = 20,
    yticklabelsvisible = true,
    yticklabelsize = 34,
    ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
)
for n in 1:Ne
        lines!(ax, μt[n], xt[n], linewidth = 2.0, color = (:blue, 0.25))
end

# Plot the bifurcation diagram
a = -3.0
b = 4.0
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
lines!(ax, μ1, x1, color = :black, linewidth = 1.5)
lines!(ax, μ2, x2, color = :black, linestyle = :dash, linewidth = 1.5)
lines!(ax, μ3, x3, color = :black, linewidth = 1.5)

# Plot the potential landscape
ax1 = Axis(fig[5:6,1],
    backgroundcolor = :transparent,
    title = L"\mu=0",
    titlesize = 34,
    xgridvisible = false,
    xlabelvisible = false,
    xticklabelsvisible = false,
    xticksvisible = false,
    ygridvisible = false,
    ylabelvisible = false,
    yticklabelsvisible = false,
    yticksvisible = false,
    limits = ((-2,5), (-2,2))
)
u=0
X = LinRange(-2.0,5.0,1000)
Y = u*X .+ X.^2 .- X.^3 .+0.2*X.^4
lines!(ax1, X, Y, color = :black, linewidth = 1.5)
for x in xt
        U = x[300000]
        scatter!(ax1, U, u*U .+ U.^2 .- U.^3 .+0.2*U.^4, color = (:blue, 0.5))
end
ax2 = Axis(fig[5:6,2],
    backgroundcolor = :transparent,
    title = L"\mu=1",
    titlesize = 34,
    xgridvisible = false,
    xlabelvisible = false,
    xticklabelsvisible = false,
    xticksvisible = false,
    ygridvisible = false,
    ylabelvisible = false,
    yticklabelsvisible = false,
    yticksvisible = false,
    limits = ((-2,4), (-1,4))
)
u=1
X = LinRange(-2.0,4.0,1000)
Y = u*X .+ X.^2 .- X.^3 .+0.2*X.^4
lines!(ax2, X, Y, color = :black, linewidth = 1.5)
for x in xt
        U = x[400000]
        scatter!(ax2, U, u*U .+ U.^2 .- U.^3 .+0.2*U.^4, color = (:blue, 0.5))
end
ax3 = Axis(fig[5:6,3],
    backgroundcolor = :transparent,
    title = L"\mu=2",
    titlesize = 34,
    xgridvisible = false,
    xlabelvisible = false,
    xticklabelsvisible = false,
    xticksvisible = false,
    ygridvisible = false,
    ylabelvisible = false,
    yticklabelsvisible = false,
    yticksvisible = false,
    limits = ((-2,4), (-1,5))
)
u=2
X = LinRange(-2.0,4.0,1000)
Y = u*X .+ X.^2 .- X.^3 .+0.2*X.^4
lines!(ax3, X, Y, color = :black, linewidth = 1.5)
for x in xt
        U = x[500000]
        scatter!(ax3, U, u*U .+ U.^2 .- U.^3 .+0.2*U.^4, color = (:blue, 0.5))
end
# Export figures
save("../results/potential/fig1.2.2.png", fig)
