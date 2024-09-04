using DynamicalSystems, DifferentialEquations
using LaTeXStrings, CairoMakie, Makie.Colors
using Statistics, ProgressMeter 

# Noise level and timescale separation
σ = 0.01
ε = 0.1
# Initial state
x0 = [3.33,-3.0]
# Temporal evolution
T = 100.00
δt = 1e-1

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
                limits = ((-3,5), (-12.5,25))
        )
        local X = LinRange(-3.0,5.0,1000)
        local Y = μ*X .+ X.^2 .- X.^3 .+0.2*X.^4
        lines!(ax, X, Y, color = :black, linewidth = 1.5)
        for x in xt
                U = x[n]
                scatter!(ax, U, μ*U .+ U.^2 .- U.^3 .+0.2*U.^4, color = (:blue, 0.5), markersize = 25)
        end
        save("../../results/precursors/B-tipping/$n.png", fig)
end
