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

                # Plot the lattice dynamical system evolution
                fig1 = Figure(; size = (1000, 1000), backgroundcolor = "#fdfff2ff")# :transparent)
                ax1 = Axis(fig1[1,1],
                          title = L"ε = %$drift,\; σ = %$noise,\; c = %$parameter",
                          backgroundcolor = "#fdfff2ff",# :transparent,
                          xticksvisible = false,
                          yticksvisible = false,
                          xticklabelsvisible = false,
                          yticklabelsvisible = false,
                          aspect = DataAspect()
                )
                heatmap!(ax1, state', colorrange = (0,1))

                # Plot the leading eigenvector evolution 
                fig2 = Figure(; size = (1000, 1000), backgroundcolor = :transparent)
                ax2 = Axis(fig2[1,1],
                          #title = L"window\; size = %$width",
                          backgroundcolor = :transparent,
                          xticksvisible = false,
                          yticksvisible = false,
                          xticklabelsvisible = false,
                          yticklabelsvisible = false,
                          aspect = DataAspect()
                )
                heatmap!(ax2, mode', colorrange = ((findmin(mode))[1],(findmax(mode))[1]))

                # Plot the leading eigenvalue evolution 
                #=
                fig3 = Figure(; size = (1000, 400))
                ax3 = Axis(fig3[1,1],
                           limits = ((0, dt*(length(eigenvalues)-width)),((findmin(eigenvalues))[1],(findmax(eigenvalues))[1])) #(0.98, 1.0))
                )
                lines!(ax3, LinRange(0.0, time, length(ews)), ews)
                =#

                # Export the figures
                save("../results/lattice/lung_ventilation/snapshots/$ctr.png", fig1)
                save("../results/lattice/lung_ventilation/modes/$ctr.png", fig2)

                # Update the time variable
                time += dt
                # Update the image counter
                ctr += 1
        end
        # Plot the leading eigenvalue evolution 
        fig3 = Figure(size = (2500, 1250), backgroundcolor = :transparent)
        ax3 = Axis(fig3[1,1],
                   backgroundcolor = :transparent,
                   aspect = 4.0,
                   xtickformat = values -> ["$(round(value*0.1,digits=2))" for value in values],
                   limits = ((0, dt*(length(eigenvalues)-width)),(0.99, 1.005))
        )
        rowsize!(fig3.layout, 1, Aspect(1, 0.25))
        resize_to_layout!(fig3)
        lines!(ax3, LinRange(0.0, time, length(eigenvalues)), eigenvalues)
        save("../results/lattice/lung_ventilation/EWS.png", fig3)
end
