"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")

# Main algorithm 
function main()
        # Parameter sweep loop
        #=@showprogress=#for (idx_μ, μ) in enumerate(μ_set[1])
                # Compute the equilibrium distribution
                N = normalise(ρ, μ)
                p(x) = N*ρ(x, μ) 

                # Compute the relevant equilibria
                equilibria = get_equilibria(f, μ, domain=[-10,10])
                xu = equilibria.unstable[1] 
                xs = equilibria.stable[2]

                # Solve the ensemble problem
                x0 = [xs, μ]
                ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

                # Ensemble loop
                for (idx_sol, solution) in enumerate(ensemble.state)
                        # Truncate and detrend the solution
                        idx_tip = find_tipping(solution, xu)
                        u = solution[1:idx_tip]
                        u = u .- x0[1]

                        # ---------------------------- # 
                        #     Parameter estimation     #
                        # ---------------------------- # 
                         
                        # E-M LLS regression
                        θ = estimate_parameters(u, dt)

                        V = fit_potential(θ, xs, U, μ)
                        error = get_error(U, V, (xs, xu), μ)
                        display(error)

                        q = fit_density(V, D)
                        error = get_error(p, q, xs, θ)
                        display(error)

                        fig = Figure()
                        ax = Axis(fig[1,1])

                        domain = LinRange(-1.5, 1.5, 100)
                        lines!(ax, domain, [p(x) for x in domain], color = :red, linewidth = 3.0)

                        domain = LinRange(0.9, 1.1, 100)
                        lines!(ax, domain, [q(x) for x in domain], color = (:blue, 0.25), linewidth = 3.0)

                        savefig("$idx_sol.png", fig)

                        println("----------------------------")
                end
        end
end

# Execute the main 
main()
