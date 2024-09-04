using DynamicalSystems, DifferentialEquations
using LaTeXStrings, CairoMakie, Makie.Colors
using Statistics, ProgressMeter 

# Noise level and timescale separation
σ = 0.00
ε = 1.1
# Initial state
x0 = [-5.0,5.0]
# Temporal evolution
T = 8.00
δt = 1e-2

# Deterministic dynamics 
function iip_det!(f, x, y, t)
        f[1] = (x[1]+x[2])^2 - 1
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
Ne = 1

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

# Plot the potential landscape
xmin = -25 
xmax = 10 
ymin = -1500 
ymax = 1500 
CairoMakie.activate!(; px_per_unit = 1)
println("Exporting the landscape evolution")
@showprogress for n in 1:length(xt[1])
        local μ = μt[1][n]
        local fig = Figure(; size = (1000, 1000), backgroundcolor = "#fdfff2ff")
        local ax = Axis(fig[1,1],
                backgroundcolor = "#fdfff2ff",
                title = L"\mu=%$(round(μ; digits = 3))",
                titlevisible = false,
                titlesize = 34,
                xgridvisible = false,
                xlabelvisible = false,
                xticklabelsvisible = false,
                xticksvisible = false,
                ygridvisible = false,
                ylabelvisible = false,
                yticklabelsvisible = false,
                yticksvisible = false,
                limits = ((xmin,xmax), (ymin,ymax))
        )
        local X = LinRange(xmin,xmax,1000)
        local Y = .-X.*(μ^2 - 1.0) .- μ.*X.^2 .-0.3*X.^3
        lines!(ax, X, Y, color = :black, linewidth = 1.5)
        for x in xt
                U = x[n]
                scatter!(ax, U, .-U.*(μ^2 - 1.0) .- μ.*U.^2 .-0.3*U.^3, color = (:green, 0.5), markersize = 25)
        end
        save("../../results/precursors/R-tipping/$n.png", fig)
end
