using DynamicalSystems
using LaTeXStrings, CairoMakie

# Define the dynamics
function iip_sn!(f, x, p, t)
        f[1] = (1/p[1])*(x[2] + x[3] + x[1]*(x[1]-1))
        f[2] = -(x[1] + (x[1])^2 + (x[1])^3 + (x[1])^4 + (x[1])^5)
        f[3] = p[2] 
        return nothing
end

# Define the set of ICs 
x0 = [[0.00, -0.50, 0.10], [1.00, -1.00, 0.10]]
# Define the parameter values
ε = 0.02
R = 1/2 + 1/4 + 1/8 + 1/16 + 1/32
r = 1.2 
p = [ε, r]
# Define the temporal parameters
T1 = 1.20
T2 = 1.20
δt = 1e-2

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
    ylabel = L"x_1(t)",
    limits = ((t[1],t[end]), (nothing,nothing))
)
lines!(ax1, t, x1, color = :green, linewidth = 1.5)
# Plot the time trajectory for y(t)
ax2 = Axis(fig1[2, 1],
    #backgroundcolor = :transparent,
    xlabel = L"t",
    ylabel = L"x_2(t)",
    limits = ((t[1],t[end]), (nothing,nothing))
)
lines!(ax2, t, x2, color = :green, linewidth = 1.5)

# Plot the trajectory in state space 
CairoMakie.activate!()
fig2 = Figure(; size = (600, 400), backgroundcolor = :transparent)
ax3 = Axis(fig2[1, 1],
    #backgroundcolor = :transparent,
    title = L"r = 1.2 > r_c = %$R",
    xlabel = L"x_1",
    ylabel = L"x_2",
    limits = ((-2.00,2.00), (-2.00,2.00))
)
lines!(ax3, x1, x2, color = :green, linewidth = 1.5)
scatter!(ax3, x0[1][1], x0[1][2], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)

# Evolve the dynamical system from the second initial condition
saddle_node = ContinuousDynamicalSystem(iip_sn!, x0[2], p)
Xt, t = trajectory(saddle_node, T2; Δt=δt)
x4 = Xt[:,1]
x5 = Xt[:,2]
x6 = Xt[:,3]

# Plot the time trajectories
lines!(ax1, t, x4, color = :orange, linewidth = 1.5)
lines!(ax2, t, x5, color = :orange, linewidth = 1.5)
# Plot the trajectory in state space
lines!(ax3, x4, x5, color = :orange, linewidth = 1.5)
scatter!(ax3, x0[2][1], x0[2][2], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)

# Plot the critical manifold and the QSE
x1_values = LinRange(-10.0,10.0,1000)
x2_values = -x3[1] .- x1_values.*(x1_values .- 1.0)
lines!(ax3, x1_values, x2_values, color = :black, linewidth = 0.75)
invariant = 0.5594.*ones(Float64, 1000)
lines!(ax3, invariant, x1_values, color = :red, linewidth = 1.0)
qse = zeros(Float64, 1000)
lines!(ax3, qse, x1_values, color = :black, linestyle = :dash, linewidth = 1.0)

# Export the figures
save("../results/CompBombTraj.png", fig1)
save("../results/CompBombState.png", fig2)

# Create animation
fig = Figure(; size = (600,400))
ax = Axis(fig[1,1],
          title = L"r = 1.2 > r_c = %$R",
          xlabel = L"x_1",
          ylabel = L"x_2",
          limits = ((-2.00,2.00), (-2.00,2.00))
)
scatter!(ax, x0[1][1], x0[1][2], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
scatter!(ax, x0[2][1], x0[2][2], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
lines!(ax, x1_values, x2_values, color = :black, linewidth = 0.75, linestyle = :dash)
lines!(ax, invariant, x1_values, color = :red, linewidth = 1.0)
lines!(ax, qse, x1_values, color = :black, linestyle = :dash, linewidth = 1.0)

frames = 1:size(x1)[1] 
framerate = 60

points1 = Observable(Point2f[(x0[1][1],x0[1][2])])
points2 = Observable(Point2f[(x0[2][1],x0[2][2])])
time = Observable(x0[1][3])
x2_values = @lift(-$time .- x1_values.*(x1_values .- 1.0))
lines!(ax, x1_values, x2_values, color = :black, linewidth = 0.75)
lambda = @lift(-$time)
scatter!(ax, 0.0, lambda, color = :black, markersize = 9)
text!(ax, 0.95, 0.95, text = @lift("λ = " .* string.(round($time + 0.0, digits = 3))))

record(fig, "../results/CompBombState.gif", frames; framerate=framerate) do frame 
       new_point1 = Point2f(x1[frame],x2[frame])
       new_point2 = Point2f(x4[frame],x5[frame])
       points1[] = push!(points1[], new_point1)
       points2[] = push!(points2[], new_point2)
       time[] = x3[frame]
       lines!(ax, points1, color = :green, linewidth = 1.5)
       lines!(ax, points2, color = :orange, linewidth = 1.5)
end
