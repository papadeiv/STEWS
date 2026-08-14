# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")

function fit_kl(u, dt; burnin = 0)
        # Rough estimate of the leading eigenvalue thriugh AR(1)
        ρ = sum(u[1:end-1].*u[2:end])/sum(abs2, u)
        λ = -log(ρ)/dt

        # Remove the transient from the input timeseries
        x = @view u[(burnin+1):end]
        Nt = length(x)

        # Compute the number of (non-overlapping) windows
        width = round(Integer, 40/(λ*dt)) + 1
        Nw = (Nt - width)÷width + 1

        # Construct the sample matrix
        sample = Matrix{Real}(undef, Nw, width)
        for idx_win in 1:Nw
                @views sample[idx_win,:] = x[((idx_win - 1)*width+1):((idx_win - 1)*width+width)]
        end
        #=
        fig = Figure()
        ax = Axis(fig[1,1], limits = (-1,1000,nothing,nothing))
        lines!(ax, LinRange(1,length(x),length(x)), x, color = :black, linewidth = 2.0)
        lines!(ax, LinRange(1,width,width), sample[1,:], color = :red, linewidth = 2.0)
        lines!(ax, LinRange(width,2*width,width), sample[2,:], color = :blue, linewidth = 2.0)
        savefig("./timeseries.png", fig)
        =#
       
        # Compute the sample covariance matrix
        K = cov(sample)

        # Solve the (symmetrised) Nystrom eigenproblem
        w  = fill(dt, width); w[1] = dt/2; w[end] = dt/2
        sw = sqrt.(w)
        S  = Symmetric((sw.*K).*sw')
        Λ = sort(eigvals(S), rev = true)

        println(" - From AR(1): $(-λ)")

        # this is where you write the bit in fit_a_sigma but in my style...
        T = (width - 1)*dt
        N = min(25, width)
        ω = [(n - 0.5)π/T for n in 1:N]
        cost(a) = (r = @. log(Λ[1:N]) + log(a^2 + ω^2); sum(abs2, r .- mean(r)))
        grid = exp.(range(log(0.1π/T), log(0.5/dt), length = 2000))
        return grid[argmin(cost.(grid))]
end

function main()
        @showprogress for (idx_μ, μ) in enumerate(μ_set)
                # Compute the relevant equilibria
                equilibria = get_equilibria(f, μ, domain=[-10,10])
                xu = equilibria.unstable[1] 
                xs = equilibria.stable[2]

                # Solve the ensemble problem
                x0 = [xs, μ]
                ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

                # Ensemble loop
                data = Matrix{Float64}(undef, convert(Integer, Ne), 4)
                for (idx_sol, solution) in enumerate(ensemble.state)
                        # Truncate and detrend the solution
                        idx_tip = find_tipping(solution, xu)
                        u = solution[1:idx_tip]
                        u = u .- x0[1]

                        # Estimate leading eigenvalue
                        println("Leading eigenvalue at μ = $μ")
                        println(" - Analytical: $(J(xs))")
                        λ = fit_kl(u, dt)
                        println(" - From KL sp: $(-λ)")
                end
        end
end

main()
