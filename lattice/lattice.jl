using CairoMakie

function neighbours(I::Tuple{Int,Int})
    neigh = Array{Int,2}(undef,4,prod(I))
    rows = I[1]
    cols = I[2]
    ctr = 0
    for ix in 1:cols
        for iy in 1:rows
            ctr += 1
            N = mod1(iy - 1, rows), ix 
            E = iy, mod1(ix + 1, cols)
            S = mod1(iy + 1, rows), ix
            W = iy, mod1(ix - 1, cols)
            neigh[1,ctr] = N[1] + I[1]*(N[2]-1)
            neigh[2,ctr] = W[1] + I[1]*(W[2]-1)
            neigh[3,ctr] = S[1] + I[1]*(S[2]-1)
            neigh[4,ctr] = E[1] + I[1]*(E[2]-1)
        end
    end
    return neigh
end

struct Lattice
        rows::Int
        cols::Int
        grid::Array{Int,2}
        connectivity::Array{Int,2}
        Lattice(rows, cols) = new(rows, cols, LinearIndices((1:rows,1:cols)), neighbours((rows,cols)))
end

function show(L::Lattice, states::Vector{Vector{Float64}}, eigenvalues::Vector{Float64}, eigenvectors::Array{Float64})
        CairoMakie.activate!()
        # Define sliding window parameters
        width = length(states) - length(eigenvalues)
        time = 0.0
        dt = 1e-2
        # Define vector storing the proposed spatial ews
        ews = Vector{Float64}() 
        # Define counter for eigenvectors 
        ctr = 1
        for j=width:(length(states)-1)
                # Get the current system's state values
                snapshot = states[j][1:end-3]
                snapshot = reshape(snapshot, (L.rows, L.cols))
                # Get the current eigenvector's state values
                eigenvec = eigenvectors[:, ctr]
                eigenvec = reshape(eigenvec, (L.rows, L.cols))
                # Get the current leading eigenvalue
                push!(ews, eigenvalues[ctr])
                # Get the model's parameters values
                parameter = round(states[j][end-2], digits=3)
                noise = states[j][end-1]
                drift = states[j][end]

                # Reshape the state and eignevector for the heatmap
                state = Array{Float64}(undef, L.rows, L.cols)
                mode = Array{Float64}(undef, L.rows, L.cols)
                for k=0:(L.rows-1)
                        state[k+1,:] = snapshot[L.rows-k,:]
                        mode[k+1,:] = eigenvec[L.rows-k,:]
                end

                fig = Figure(; size = (600, 400))
                # Plot the lattice dynamical system evolution
                ax1 = Axis(fig[1, 1],
                          title = L"ε = %$drift,\; σ = %$noise,\; k = %$parameter",
                          xticksvisible = false,
                          yticksvisible = false,
                          xticklabelsvisible = false,
                          yticklabelsvisible = false,
                          aspect = DataAspect()
                )
                
                heatmap!(ax1, state', colorrange = (0,1))
                # Plot the leading eigenvector evolution 
                ax2 = Axis(fig[1, 2],
                          title = L"window\; size = %$width",
                          xticksvisible = false,
                          yticksvisible = false,
                          xticklabelsvisible = false,
                          yticklabelsvisible = false,
                          aspect = DataAspect()
                )
                #=
                Colorbar(fig[:, end+1], 
                         ticks = [0,1]
                )
                =#
                heatmap!(ax2, mode', colorrange = ((findmin(mode))[1],(findmax(mode))[1]))
                # Plot the leading eigenvalue evolution 
                ax3 = Axis(fig[2, 1:2],
                           limits = ((0, dt*(length(eigenvalues)-width)), (0.98, 1.0))
                )
                lines!(ax3, LinRange(0.0, time, length(ews)), ews)

                # Export the figure
                save("../results/lattice/lung_ventilation/$ctr.png", fig)

                # Update the time variable
                time += dt
                # Update the image counter
                ctr += 1
        end
end
