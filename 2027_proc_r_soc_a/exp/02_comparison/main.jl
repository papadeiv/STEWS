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
        # Materialise the simulation parameter grid for parallelization
        grid = vec(collect(Iterators.product(dt, Nt)))

        # Pre-create the output directories to avoid conflicts among threads
        for idx_μ in eachindex(μ_set), (stepsize, _) in grid
                mkpath("../../res/data/$idx_μ/$stepsize")
        end

        # Simulation parameters loop
        printstyled("Generating the samples using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
        progressbar = Progress(length(grid))
        Threads.@threads :dynamic for idx_sim in eachindex(grid) 
                # Specify the simulation setup
                stepsize, steps = grid[idx_sim]

                # Solve the ensemble problems and export the results
                tipped = generate_samples(stepsize, steps)
                next!(progressbar)
        end

        # Perform the statistical analysis of the results 
        analysis()
        printstyled("Simulation completed!\n"; bold=true, color=:green)
end

# Execute the main 
main()
