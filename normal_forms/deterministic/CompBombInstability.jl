using DynamicalSystems
using LaTeXStrings, CairoMakie

# Define the dynamics
function iip_sn!(f, x, p, t)
        f[1] = (1/p[1])*(x[2] + x[3] + x[1]*(x[1]-1))
        f[2] = -(x[1] + x[1]^2 + x[1]^3)
        f[3] = p[2] 
        return nothing
end

# Define the set of ICs 
x0 = [[0.00, -0.50, -1.00], [0.00, 0.50, 1.00]]
# Define the parameter values
epsilon = 0.02
r = 0.50
p = [epsilon, r]
# Define the temporal parameters
T1 = 20.00
T2 = 20.00
δt = 1e-3

# Evolve the dynamical system from the first initial condition
saddle_node = ContinuousDynamicalSystem(iip_sn!, x0[1], p)
Xt, t = trajectory(saddle_node, T1; Δt=δt)
x1 = Xt[:,1]
x2 = Xt[:,2]
x3 = Xt[:,3]

# Plot the time trajectory for x(t)
CairoMakie.activate!()
fig1 = Figure(; size = (600, 400), backgroundcolor = :transparent)
ax1 = Axis(fig1[1, 1],
    #backgroundcolor = :transparent,
    xlabel = L"t",
    ylabel = L"x(t)",
    limits = ((t[1],t[end]), (nothing,nothing))
)
lines!(ax1, t, x, color = :red, linewidth = 1.5)
# Plot the time trajectory for y(t)
ax2 = Axis(fig1[2, 1],
    #backgroundcolor = :transparent,
    xlabel = L"t",
    ylabel = L"y(t)",
    limits = ((t[1],t[end]), (nothing,nothing))
)
lines!(ax2, t, y, color = :blue, linewidth = 1.5)

# Plot the trajectory in state space 
CairoMakie.activate!()
fig2 = Figure(; size = (600, 400), backgroundcolor = :transparent)
ax3 = Axis(fig2[1, 1],
    #backgroundcolor = :transparent,
    title = L"r = %$r",
    xlabel = L"x",
    ylabel = L"y",
    limits = ((-2.00,2.00), (-2.00,2.00))
)
lines!(ax3, x, y, color = :green, linewidth = 1.5)
scatter!(ax3, x0[1][1], x0[1][2], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)

# Evolve the dynamical system from the second initial condition
saddle_node = ContinuousDynamicalSystem(iip_sn!, x0[2], p)
Xt, t = trajectory(saddle_node, T2; Δt=δt)
x1 = Xt[:,1]
x2 = Xt[:,2]
x3 = Xt[:,3]

# Plot the time trajectories
lines!(ax1, t, x, color = :red, linewidth = 1.5, linestyle = :dash)
lines!(ax2, t, y, color = :blue, linewidth = 1.5, linestyle = :dash)
# Plot the trajectory in state space
lines!(ax3, x, y, color = :yellow, linewidth = 1.5)
scatter!(ax3, x0[2][1], x0[2][2], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)

# Plot the critical manifold
y_values = LinRange(-10.0,10.0,10)
attracting = sqrt(mu-r) .- y_values
repelling = -sqrt(mu-r) .- y_values
lines!(ax3, y_values, attracting, color = :blue, linewidth = 0.75)
lines!(ax3, y_values, repelling, color = :red, linewidth = 0.75)

#Plot the isoclines
unstable = sqrt(mu) .- y_values
stable = -sqrt(mu) .- y_values
lines!(ax3, y_values, stable, color = :black, linewidth = 0.75, linestyle = :dash)
lines!(ax3, y_values, unstable, color = :black, linewidth = 0.75, linestyle = :dash)

# Export the figures
save("../results/CompBombTraj.png", fig1)
save("../results/CompBombState.png", fig2)
