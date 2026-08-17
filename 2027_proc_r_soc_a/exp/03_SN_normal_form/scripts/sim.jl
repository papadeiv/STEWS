"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ_set = collect(range(-0.9, stop=-0.2, step=0.1)) # Bifurcation parameter value
ε = 0.0                                           # Timescale separation
σ = 0.010                                         # Noise level (additive)
D = (σ^2)/2.0                                     # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ - x^2                                # Drift
Λ(t) = ε                                          # Shift/Ramp
η(x) = σ                                          # Diffusion

# Simulation parameters
dt = 1e-1                                         # Timestep
Nt = 1e3                                          # Total number of steps
Ne = 1e2                                          # Number of particles in the ensemble 

function estimate_parameters_2024(x::AbstractVector{<:Real}, dt=1.0; nbins=10, K::Int=0, iters=8)
    edges = range(minimum(x), maximum(x), nbins+1); ctr = (edges[1:end-1].+edges[2:end])./2
    bin(v) = clamp(searchsortedlast(edges,v), 1, nbins)
    n=zeros(nbins); Sy=zeros(nbins); Syy=zeros(nbins)
    T=zeros(nbins,K); Ry=zeros(nbins,K); U=zeros(nbins,K,K)
    for i in K+1:length(x)-1                               # sufficient statistics
        b=bin(x[i]); yv=x[i+1]-x[i]; Δ=[x[i]-x[i-k] for k in 1:K]
        n[b]+=1; Sy[b]+=yv; Syy[b]+=yv^2
        T[b,:].+=Δ; Ry[b,:].+=yv.*Δ; U[b,:,:].+=Δ*Δ'
    end
    act=findall(>(0),n); na=length(act); p=na+K
    v=[n[b]>0 ? Syy[b]/n[b] : NaN for b in 1:nbins]; D1=fill(NaN,nbins); κ=zeros(K)
    for _ in 1:iters                                       # GLS + variance update
        w=[n[b]>0 ? 1/v[b] : 0.0 for b in 1:nbins]; G=zeros(p,p); h=zeros(p)
        for (a,b) in enumerate(act)
            G[a,a]=dt^2*w[b]*n[b]; h[a]=dt*w[b]*Sy[b]
            for k in 1:K; G[a,na+k]=G[na+k,a]=dt^2*w[b]*T[b,k]; end
        end
        for k in 1:K, l in 1:K; G[na+k,na+l]=dt^2*sum(w[b]*U[b,k,l] for b in act); end
        for k in 1:K; h[na+k]=dt*sum(w[b]*Ry[b,k] for b in act); end
        θ=G\h; κ=θ[na+1:end]
        for (a,b) in enumerate(act)
            d=θ[a]; D1[b]=d
            v[b]=(Syy[b]-2dt*(d*Sy[b]+κ⋅Ry[b,:])+dt^2*(d^2*n[b]+2d*(κ⋅T[b,:])+κ⋅(U[b,:,:]*κ)))/n[b]
        end
    end
    i = findfirst(i->D1[i]>0>D1[i+1], 1:nbins-1)           # stable zero
    i===nothing && return (ctr, D1, NaN)
    xs = ctr[i] - D1[i]*(ctr[i+1]-ctr[i])/(D1[i+1]-D1[i])
    m = max(1,i-1):min(nbins,i+2); w=n[m]; c=ctr[m].-xs; r=D1[m]  # count-weighted slope
    A = [sum(w) sum(w.*c); sum(w.*c) sum(w.*c.^2)] \ [sum(w.*r); sum(w.*c.*r)]
    (grid=ctr, D1=D1, α=-A[2])
end

function estimate_parameters_2009(x::AbstractVector{<:Real}, dt; ngrid=50, alpha_bw=nothing)
    κ(u) = u^2 < 5 ? 1 - u^2/5 : 0.0
    xj, y, τ = x[1:end-1], diff(x), dt
    m = min(length(xj), 3000)
    h0 = std(x)*length(x)^(-0.2); h, best = h0, Inf
    for hc in range(0.3h0, 3h0, 12)                     # CV bandwidth
        s = 0.0
        for k in 1:m
            n = d = 0.0
            for j in 1:m
                j == k && continue
                w = κ((xj[k]-xj[j])/hc); n += w*y[j]; d += w
            end
            d > 0 && (s += (y[k]/τ - n/(d*τ))^2)
        end
        s < best && ((best, h) = (s, hc))
    end
    grid = range(quantile(x,0.02), quantile(x,0.98), ngrid)
    function nw(x0, b)                                  # Nadaraya–Watson drift at x0
        s = w = 0.0
        for j in eachindex(xj)
            k = κ((x0-xj[j])/b); w += k; s += k*y[j]/τ
        end
        w > 0 ? s/w : NaN
    end
    D1 = nw.(grid, h)
    i = findfirst(i -> D1[i] > 0 > D1[i+1], 1:ngrid-1)  # stable zero
    i === nothing && return (grid, D1, NaN)
    xs = grid[i] - D1[i]*step(grid)/(D1[i+1]-D1[i])
    hl = alpha_bw === nothing ? h : alpha_bw            # local-linear slope at xs
    c = xj .- xs; w = κ.(c ./ hl); r = y ./ τ
    A = [sum(w) sum(w.*c); sum(w.*c) sum(w.*c.^2)] \ [sum(w.*r); sum(w.*c.*r)]
    return (grid=grid, D1=D1, α=-A[2])
end
