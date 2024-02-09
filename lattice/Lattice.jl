using CairoMakie

function neighbours(I::Tuple{Int,Int})
    neigh = Array{Int,2}(undef,4,prod(I))
    Ny = I[1] # Number of rows
    Nx = I[2] # Number of columns
    ctr = 0
    for ix in 1:Nx
        for iy in 1:Ny        
            ctr += 1
            N = mod1(iy - 1, Ny), ix 
            E = iy, mod1(ix + 1, Nx)
            S = mod1(iy + 1, Ny), ix
            W = iy, mod1(ix - 1, Nx)
            neigh[1,ctr] = N[1] + I[1]*(N[2]-1)
            neigh[2,ctr] = W[1] + I[1]*(W[2]-1)
            neigh[3,ctr] = S[1] + I[1]*(S[2]-1)
            neigh[4,ctr] = E[1] + I[1]*(E[2]-1)
        end
    end
    return neigh
end

struct Lattice
        Nx::Int
        Ny::Int
        grid::Array{Int,2}
        connectivity::Array{Int,2}
        Lattice(Nx, Ny) = new(Nx, Ny, LinearIndices((1:Nx,1:Ny)), neighbours((Ny,Nx)))
end

function show(L::Lattice, states::Vector{Vector{Float64}})
        CairoMakie.activate!()
        for j=1:length(states)
                state = states[j]
                state = reshape(state, L.Nx, L.Ny)

                fig = Figure(; size = (600, 400))
                ax = Axis(fig[1, 1],
                          xticksvisible = false,
                          yticksvisible = false,
                          xticklabelsvisible = false,
                          yticklabelsvisible = false 
                )
                Colorbar(fig[:, end+1], 
                         ticks = [0,1]
                )
                heatmap!(ax, 1:L.Nx, 1:L.Ny, state, colorrange = (0,0.1))#highclip = :green, lowclip = :red, colorrange = (0.0,1.0))

                save("../results/lattice/lung_ventilation/$j.png", fig)
        end
end
