using DynamicalSystems, DifferentialEquations
using LaTeXStrings, CairoMakie, Makie.Colors
using Statistics, ProgressMeter 

printstyled("FOLD NORMAL FORM\n"; bold=true, underline=true, color=:light_blue)

# Define the parameters of the process
σ = 0.1
ε = 0.02

# Define initial states
x0 = [1.0,-1.0]

# Define temporal evolution quantities
T = 50.00
δt = 1e-3

# Define the transcritical form
function iip_det!(f, x, y, t)
        f[1] = -x[2] - (x[1])^2
        f[2] = ε
        return nothing
end
function iip_stoc!(f, x, y, t)
        f[1] = +σ
        f[2] = 0.0
        return nothing
end
normal_form = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T))

# Define ensemble quantities
Ne = 1000
N_esc = 0
N_uns = 0
critical_y = Float64[]

# Define parameter values
Ny = 100
y_values = LinRange(x0[2],0,Ny)

# Define the 'trajectory is escaped' termination condition
function condition_1(x, t, integrator)
        x[1] < -sqrt(abs(x[2]))
end
function condition_2(x, t, integrator)
        x[1] < 5*x[2] - 0.5
end
function condition_3(x, t, integrator)
        x[1] < -1.0 
end
function affect!(integrator)
        terminate!(integrator)
end
cb_1 = DiscreteCallback(condition_1, affect!, save_positions=(false,false))
cb_2 = DiscreteCallback(condition_2, affect!, save_positions=(false,false))
cb_3 = DiscreteCallback(condition_3, affect!, save_positions=(false,false))

# Define arrays to store the ensemble sample paths for escaped and unescaped trajectories 
xt = fill(Float64[], 0) 
yt = fill(Float64[], 0)
xt_esc = fill(Float64[], 0)
yt_esc = fill(Float64[], 0)

println("Simulating the ensemble sample paths")
# Loop over the ensemble sample paths
@showprogress for j in 1:Ne
        local sol = solve(normal_form, EM(), dt=δt, callback=cb_1, verbose=false)
        if sol[1,end] < -sqrt(abs(sol[2,end])) # It means that this trajectory has escaped the manifold 
                # Store the times series up to the escape
                push!(xt_esc, sol[1,:])
                push!(yt_esc, sol[2,:])
                # Store the value of the slow variable at the unstable critical submanifold
                push!(critical_y, sol[2,end])
                # Update the counter of escaped trajectories
                global N_esc = N_esc + 1
        else # It means that the trajectory stayed bounded for the whole time
                push!(xt, sol[1,:])
                push!(yt, sol[2,:])
                # Update the counter of unescaped trajectories
                global N_uns = N_uns + 1
        end
end
# Plot the ensemble's sample paths
CairoMakie.activate!(; px_per_unit = 3)
fig1 = Figure(; size = (1200, 400), backgroundcolor = :transparent)
ax1 = Axis(fig1[1,1:2],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((x0[2],y_values[end]), (-1.25,+1.25)),
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
    xticks = [-1,0],
    xticksvisible = true,
    xticksize = 6,
    xticklabelsvisible = true,
    xticklabelsize = 18,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"x",
    ylabelvisible = true,
    ylabelsize = 20,
    ylabelcolor = :black,
    ylabelpadding = -25.0,
    yticks = [-1,1],
    yticksvisible = true,
    yticksize = 6,
    yticklabelsvisible = true,
    yticklabelsize = 18,
    ytickformat = "{:.1f}",
    yscale = identity,
    yaxisposition = :left,

)
ax2 = Axis(fig1[2,1:2],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((x0[2],y_values[end]), (-1.25,+1.25)),
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
    xticks = [-1,0],
    xticksvisible = true,
    xticksize = 6,
    xticklabelsvisible = true,
    xticklabelsize = 18,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"x",
    ylabelvisible = true,
    ylabelsize = 20,
    ylabelcolor = :black,
    ylabelpadding = -25.0,
    yticks = [-1,1],
    yticksvisible = true,
    yticksize = 6,
    yticklabelsvisible = true,
    yticklabelsize = 18,
    ytickformat = "{:.1f}",
    yscale = identity,
    yaxisposition = :left,
)
# Unescaped trajectories
for n in 1:N_uns # == size(xt,1)
        lines!(ax1, yt[n], xt[n], linewidth = 0.5, color = (:gray,0.15))
end
# Escaped trajectories
for n in 1:N_esc # == size(xt_esc,1)
        lines!(ax2, yt_esc[n], xt_esc[n], linewidth = 0.5, color = (:red,0.25))
end
# Plot the (deterministic) critical manifold
stable = sqrt.(-y_values)
unstable = -sqrt.(-y_values)
lines!(ax2, y_values, stable, color = :black, linewidth = 1.5)
lines!(ax2, y_values, unstable, color = :black, linewidth = 1.5, linestyle = :dash)
lines!(ax1, y_values, unstable, color = :black, linewidth = 1.5, linestyle = :dash)
# Plot percentage of escaped trajectories at different y values
ax3 = Axis(fig1[1:2,3],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((-1,0), (-0.5,100)),
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
    xticks = [-1,0],
    xticksvisible = true,
    xticksize = 6,
    xticklabelsvisible = true,
    xticklabelsize = 20,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"\text{N-tipping (%)}",
    ylabelvisible = true,
    ylabelsize = 18,
    ylabelcolor = :black,
    ylabelpadding = -25.0,
    yticks = [0,100],
    yticksvisible = true,
    yticksize = 6,
    yticklabelsvisible = true,
    yticklabelsize = 20,
    ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
)
pct_esc = Float64[]
push!(pct_esc, 0.0)
counter = 0
println("Computing the parameter distribution of the escaped trajectories")
# Loop over the discretised parameter values
@showprogress for n in 2:Ny
        for m in 1:N_esc
                w = yt_esc[m]
                if w[end] > y_values[n-1] && w[end] < y_values[n] 
                        global counter = counter + 1
                end
        end
        push!(pct_esc, (counter/Ne)*100)
end
scatter!(ax3, y_values, pct_esc, color = :red, strokecolor = :black, strokewidth = 1.5, markersize = 7)
text!(-0.75, 70, text="Tot. esc. = " .*string.(round((N_esc/Ne)*100, digits=2)) .* "%", align = (:center, :center))
# Plot the ensemble moments of the distributions for different parameter values
M = Float64[]
V = Float64[]
println("Computing the variance of the (unescaped) ensemble paths")
# Loop over the parameter values (timesteps)
@showprogress for n in 1:size(xt[end],1) 
        # Compute the ensemble mean
        MEAN = 0.0
        # Loop over the unescaped ensemble paths all at the n-th timestep
        for m in 1:N_uns # == size(xt,1)
                u = xt[m]
                MEAN = MEAN + u[n] 
        end
        # Loop over the escaped ensemble paths all at the n-th timestep
        N_esc_crit = 0
        for m in 1:N_esc # == size(xt_esc,1)
                u = xt_esc[m]
                w = yt_esc[m]
                if n <= size(w,1)
                        if w[n] < critical_y[m] 
                                MEAN = MEAN + u[n]
                                N_esc_crit = N_esc_crit + 1
                        end
                end
        end
        MEAN = MEAN/(N_uns + N_esc_crit)
        # Compute the ensemble variance
        VAR = 0.0
        # Loop over the unescaped ensemble paths all at the n-th timestep
        for m in 1:N_uns
                u = xt[m]
                VAR = VAR + (u[n] - MEAN)^(2.0) 
        end
        # Loop over the escaped ensemble paths all at the n-th timestep
        N_esc_crit = 0
        for m in 1:N_esc # == size(xt_esc,1)
                u = xt_esc[m]
                w = yt_esc[m]
                if n <= size(w,1)
                        if w[n] < critical_y[m] 
                                VAR = VAR + (u[n] - MEAN)^(2.0)
                                N_esc_crit = N_esc_crit + 1
                        end
                end
        end
        VAR = VAR/(N_uns + N_esc_crit)
        # Store the moments of the ensemble at the n-th timestep
        push!(M, MEAN)
        push!(V, VAR)
end
lines!(ax1, yt[end], M, color = :red, linewidth = 3.5)
fig2 = Figure(; size = (1200, 700), backgroundcolor = "#eeeeeeff")
ax4 = Axis(fig2[1,1],
    # Background
    backgroundcolor = :white,
    spinewidth = 5.0,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((-1,0), nothing),
    # Title
    #title = L"r = r_c = %$R",
    titlevisible = false,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"\mathbf{\mu}",
    xlabelvisible = false,
    xlabelsize = 50,
    xlabelcolor = :black,
    xlabelpadding = -10.0,
    xticks = [-1,0],
    xticksvisible = true,
    xticksize = 20,
    xtickwidth = 5.0,
    ytickwidth = 5.0,
    xticklabelsvisible = false,
    xticklabelsize = 50,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"\textbf{variance}",
    ylabelvisible = true,
    ylabelsize = 50,
    ylabelcolor = :black,
    ylabelpadding = -60.0,
    yticks = [0,0.0175],
    yticksvisible = true,
    yticksize = 20,
    yticklabelsvisible = true,
    yticklabelsize = 50,
    ytickformat = "{:.2f}",
    yscale = identity,
    yaxisposition = :left,
)
# Plot the theoretical variance of the steady-state Fokker-Planck distribution of the normal form
theoretical_y = LinRange(x0[2],0,200)
theoretical_V = [0.0025063,0.00251264,0.00251902,0.00252546,0.00253194,0.00253848,0.00254506,0.0025517,0.00255839,0.00256513,0.00257193,0.00257878,0.00258569,0.00259265,0.00259967,0.00260675,0.00261388,0.00262107,0.00262833,0.00263564,0.00264302,0.00265046,0.00265796,0.00266553,0.00267316,0.00268085,0.00268862,0.00269645,0.00270435,0.00271233,0.00272037,0.00272848,0.00273667,0.00274494,0.00275327,0.00276169,0.00277018,0.00277875,0.00278741,0.00279614,0.00280496,0.00281386,0.00282285,0.00283192,0.00284108,0.00285034,0.00285968,0.00286912,0.00287865,0.00288827,0.002898,0.00290782,0.00291775,0.00292777,0.00293791,0.00294815,0.00295849,0.00296895,0.00297952,0.0029902,0.003001,0.00301192,0.00302296,0.00303412,0.0030454,0.00305682,0.00306836,0.00308004,0.00309185,0.0031038,0.00311589,0.00312812,0.0031405,0.00315302,0.0031657,0.00317854,0.00319153,0.00320468,0.003218,0.00323149,0.00324515,0.00325898,0.003273,0.00328719,0.00330158,0.00331615,0.00333093,0.0033459,0.00336108,0.00337647,0.00339207,0.00340789,0.00342394,0.00344022,0.00345673,0.00347348,0.00349049,0.00350775,0.00352526,0.00354305,0.00356111,0.00357945,0.00359807,0.003617,0.00363623,0.00365577,0.00367563,0.00369582,0.00371636,0.00373724,0.00375848,0.00378009,0.00380208,0.00382446,0.00384725,0.00387045,0.00389408,0.00391815,0.00394268,0.00396769,0.00399318,0.00401917,0.00404568,0.00407273,0.00410034,0.00412852,0.0041573,0.0041867,0.00421674,0.00424744,0.00427884,0.00431095,0.0043438,0.00437743,0.00441186,0.00444713,0.00448328,0.00452033,0.00455834,0.00459734,0.00463737,0.00467849,0.00472074,0.00476418,0.00480886,0.00485485,0.00490222,0.00495104,0.00500137,0.00505331,0.00510694,0.00516236,0.00521968,0.005279,0.00534045,0.00540416,0.00547028,0.00553897,0.00561041,0.00568479,0.00576232,0.00584325,0.00592783,0.00601636,0.00610918,0.00620666,0.00630923,0.00641736,0.00653163,0.00665267,0.00678124,0.00691823,0.00706471,0.00722197,0.00739157,0.00757547,0.00777613,0.00799662,0.00824088,0.00851387,0.00882185,0.00917253,0.00957509,0.0100399,0.0105779,0.0111987,0.0119078,0.0127029,0.0135685,0.0144707,0.0153535,0.0161394,0.0167353,0.0170465,0.0169933,0.0165281,0.0156426,0.0143638,0.0127298,0.0107238]
lines!(ax4, theoretical_y, theoretical_V, color = :red, linewidth = 5)
lines!(ax4, yt[end], V, color = :black, linewidth = 3.5)
#text!(-0.85, 0.0155, text = L"\text{(a)}", align = [:right, :bottom], color = :black, fontsize = 50)
# Export figures
save("../fig/fig2.6.1.png", fig1)
save("../fig/fig2.5.1.png", fig2)

# Plot an additional figure for the presentation
fig3 = Figure(; size = (1200, 700), backgroundcolor = "#eeeeeeff")
ax = Axis(fig3[1,1],
    # Background
    backgroundcolor = :white,
    spinewidth = 5.0,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((x0[2],y_values[end]), (-1.25,+1.25)),
    # Title
    #title = L"r = r_c = %$R",
    titlevisible = false,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"\mathbf{\mu}",
    xlabelvisible = false,
    xlabelsize = 50,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    xticks = [-1,0],
    xticksvisible = true,
    xticksize = 20,
    xticklabelsvisible = false,
    xticklabelsize = 20,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"\mathbf{x}",
    ylabelvisible = true,
    ylabelsize = 50,
    ylabelcolor = :black,
    ylabelpadding = -60.0,
    yticks = [-1,1],
    yticksvisible = true,
    yticksize = 20,
    yticklabelsvisible = true,
    yticklabelsize = 50,
    ytickformat = "{:.1f}",
    yscale = identity,
    xtickwidth = 5.0,
    ytickwidth = 5.0,
    yaxisposition = :left,
)
for n in 1:N_uns # == size(xt,1)
        lines!(ax, yt[n], xt[n], linewidth = 0.5, color = (:gray,0.85))
end
# Escaped trajectories
for n in 1:N_esc # == size(xt_esc,1)
        lines!(ax, yt_esc[n], xt_esc[n], linewidth = 0.5, color = (:red,0.15))
end
# Plot the (deterministic) critical manifold
stable = sqrt.(-y_values)
unstable = -sqrt.(-y_values)
lines!(ax, y_values, stable, color = :black, linewidth = 5)
lines!(ax, y_values, unstable, color = :black, linewidth = 5, linestyle = :dash)
lines!(ax, y_values, unstable, color = :black, linewidth = 5, linestyle = :dash)
save("../fig/slide6.1.png", fig3)
