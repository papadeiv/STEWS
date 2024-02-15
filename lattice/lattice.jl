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

function show(L::Lattice, states::Vector{Vector{Float64}})
        CairoMakie.activate!()
        for j=1:length(states)
                snapshot = states[j][1:end-3]
                snapshot = reshape(snapshot, (L.rows, L.cols))

                parameter = round(states[j][end-2], digits=3)
                noise = states[j][end-1]
                drift = states[j][end]

                state = Array{Float64}(undef, L.rows, L.cols)
                for k=0:(L.rows-1)
                        state[k+1,:] = snapshot[L.rows-k,:]
                end

                fig = Figure(; size = (600, 400))
                ax = Axis(fig[1, 1],
                          title = L"ε = %$drift,\; σ = %$noise,\; k = %$parameter",
                          xticksvisible = false,
                          yticksvisible = false,
                          xticklabelsvisible = false,
                          yticklabelsvisible = false,
                          aspect = DataAspect()
                )
                Colorbar(fig[:, end+1], 
                         ticks = [0,1]
                )
                heatmap!(ax, state', colorrange = (0,1))

                save("../results/lattice/lung_ventilation/$j.png", fig)
        end
end
