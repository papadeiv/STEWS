using DifferentialEquations, NetworkDynamics, Graphs

# https://pik-icone.github.io/NetworkDynamics.jl/dev/getting_started_with_network_dynamics/

function diffusionedge!(e, v_s, v_d, p, t)
    # usually e, v_s, v_d are arrays, hence we use the broadcasting operator .
    e .= v_s - v_d
    nothing
end

function diffusionvertex!(dv, v, edges, p, t)
    # usually v, edges are arrays, hence we use the broadcasting operator .
    dv .= 0.0
    for e in edges
        dv .+= e
    end
    nothing
end

N = 20 # number of nodes
k = 4  # average degree
g = barabasi_albert(N, k) # a little more exciting than a bare random graph 

nd_diffusion_vertex = ODEVertex(; f=diffusionvertex!, dim=1)
nd_diffusion_edge = StaticEdge(; f=diffusionedge!, dim=1)
nd = network_dynamics(nd_diffusion_vertex, nd_diffusion_edge, g)

x0 = randn(N) # random initial conditions
ode_prob = ODEProblem(nd, x0, (0.0, 4.0))
sol = solve(ode_prob, Tsit5());
