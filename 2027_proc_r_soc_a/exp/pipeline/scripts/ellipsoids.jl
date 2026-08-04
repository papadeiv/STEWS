# ellipsoids.jl
#
# Visualise the 3-D covariance ellipsoid of a cloud of least-squares solutions
# θ = (θ₁, θ₂, θ₃) and its projections onto the three coordinate planes, for a
# family of potentials parameterised by μ ∈ [0, 1].
#
# Requires: CairoMakie ≥ 0.11 (Makie ≥ 0.20).  On older Makie replace
# `size = ...` by `resolution = ...` and `shading = NoShading` by `shading = false`.
#
#   julia> include("ellipsoids.jl")
#   julia> Xs, mus = demo_data();               # ← replace with your own data
#   julia> fig = plot_ellipsoid_family(Xs, mus)
#   julia> save("ellipsoids.png", fig)
#
# `Xs` is a vector of 100×3 matrices (one per value of μ), `mus` the matching
# vector of μ values.  Row i of Xs[k] is the solution θ of the LLS problem for
# realisation i at μ = mus[k].

using CairoMakie
using LinearAlgebra
using Statistics
using Random

CairoMakie.activate!(type = "png", px_per_unit = 2)

"""
    ellipse2d(c, Σ, dims, r; kind = :shadow, n = 181)

Boundary of the 2-D ellipse in the coordinate plane `dims` (e.g. `(1, 3)`).

* `kind = :shadow`  — the *orthogonal projection* (silhouette) of the 3-D
  ellipsoid.  It is exactly the confidence ellipse of the marginal 2×2 block
  `Σ[dims, dims]`, at the same radius `r`.  This is what you almost always want.
* `kind = :section` — the *central cross-section*, i.e. the ellipse you get by
  cutting with the plane through `c`.  It uses `inv(inv(Σ)[dims, dims])`, the
  conditional covariance, and is strictly smaller than the shadow.
"""
function ellipse2d(c, Σ, dims::NTuple{2,Int}, r; kind = :shadow, n = 181)
    idx = collect(dims)
    S = kind === :shadow  ? Σ[idx, idx] :
        kind === :section ? inv(inv(Symmetric(Matrix(Σ)))[idx, idx]) :
        throw(ArgumentError("kind must be :shadow or :section"))
    L = sqrt_cov(S)
    t  = range(0, 2π; length = n)
    xs = Vector{Float64}(undef, n); ys = similar(xs)
    for (i, τ) in enumerate(t)
        p = c[idx] .+ r .* (L * [cos(τ), sin(τ)])
        xs[i], ys[i] = p
    end
    return xs, ys
end

# ----------------------------------------------------------------------------
# 2.  The figure
# ----------------------------------------------------------------------------

"""
    plot_ellipsoid_family(Xs, mus; kwargs...)

`Xs::Vector{<:AbstractMatrix}` — one `100×3` matrix per value of μ.
`mus::AbstractVector`          — the μ values, used for the colour code.

Keywords
* `conf = 0.95`         confidence level of the ellipsoids
* `style = :surface`    `:surface` (translucent, few μ) or `:wire` (many μ)
* `alpha = 0.22`        opacity of the surfaces
* `kind = :shadow`      `:shadow` or `:section` for the 2-D panels
* `points = :ends`      scatter the raw solutions: `:none`, `:ends`, `:all`
* `shadows = true`      also draw the three projections on the walls of the box
* `colormap = :viridis`
* `labels = ("θ₁","θ₂","θ₃")`
"""
function plot_ellipsoid_family(Xs::AbstractVector{<:AbstractMatrix},
                               mus::AbstractVector;
                               conf     = 0.95,
                               style    = :surface,
                               alpha    = 0.22,
                               kind     = :shadow,
                               points   = :ends,
                               shadows  = true,
                               colormap = :viridis,
                               labels   = ("θ₁", "θ₂", "θ₃"))

    length(Xs) == length(mus) || throw(ArgumentError("Xs and mus must match"))
    ells = [fit_ellipsoid(X; conf) for X in Xs]          # (c, Σ, r) per μ

    # colour code -------------------------------------------------------------
    grad = cgrad(colormap)
    μlo, μhi = extrema(mus)
    colof(μ) = grad[μhi > μlo ? (μ - μlo) / (μhi - μlo) : 0.5]

    # global bounding box: the support function of an ellipsoid along eₖ is
    # cₖ ± r·√Σₖₖ, so this box is tight.
    lo = fill(Inf, 3); hi = fill(-Inf, 3)
    for (c, Σ, r) in ells, k in 1:3
        lo[k] = min(lo[k], c[k] - r*sqrt(Σ[k, k]))
        hi[k] = max(hi[k], c[k] + r*sqrt(Σ[k, k]))
    end
    pad = 0.06 .* (hi .- lo); lo .-= pad; hi .+= pad

    fig = Figure(size = (1150, 950))

    ax3 = Axis3(fig[1:3, 1];
                xlabel = labels[1], ylabel = labels[2], zlabel = labels[3],
                title  = "$(round(Int, 100conf))% ellipsoids",
                azimuth = 1.275π, elevation = π/8,      # wall choice below assumes this
                limits = (lo[1], hi[1], lo[2], hi[2], lo[3], hi[3]))

    planes = ((1, 2), (1, 3), (2, 3))
    axs2 = [Axis(fig[i, 2];
                 xlabel = labels[p[1]], ylabel = labels[p[2]],
                 title  = "($(labels[p[1]]), $(labels[p[2]])) plane")
            for (i, p) in enumerate(planes)]

    # which wall each shadow is glued to, for azimuth = 1.275π (camera in the
    # −x, −y octant ⇒ far walls are +x, +y, floor is −z)
    walls = ((1, 2) => (3, lo[3]), (1, 3) => (2, hi[2]), (2, 3) => (1, hi[1]))

    for (k, ((c, Σ, r), μ)) in enumerate(zip(ells, mus))
        col  = colof(μ)
        show = points === :all || (points === :ends && (k == 1 || k == length(ells)))

        # ---- 3-D ellipsoid --------------------------------------------------
        if style === :surface
            Xg, Yg, Zg = ellipsoid_surface(c, Σ, r)
            surface!(ax3, Xg, Yg, Zg;
                     color = fill(RGBAf(col.r, col.g, col.b, alpha), size(Xg)),
                     shading = NoShading, transparency = true)
        else
            lines!(ax3, ellipsoid_wireframe(c, Σ, r);
                   color = (col, 0.7), linewidth = 0.9)
        end

        # ---- shadows on the walls of the 3-D box ----------------------------
        if shadows
            for (p, (kfix, val)) in walls
                xs, ys = ellipse2d(c, Σ, p, r; kind = :shadow)
                P = map(eachindex(xs)) do i
                    q = zeros(3); q[p[1]] = xs[i]; q[p[2]] = ys[i]; q[kfix] = val
                    Point3f(q)
                end
                lines!(ax3, P; color = (col, 0.35), linewidth = 1)
            end
        end

        # ---- 2-D projections ------------------------------------------------
        for (ax, p) in zip(axs2, planes)
            xs, ys = ellipse2d(c, Σ, p, r; kind)
            poly!(ax, Point2f.(xs, ys);
                  color = (col, 0.12), strokecolor = col, strokewidth = 1.6)
            show && scatter!(ax, Xs[k][:, p[1]], Xs[k][:, p[2]];
                             color = (col, 0.45), markersize = 4)
            scatter!(ax, [c[p[1]]], [c[p[2]]];
                     color = col, markersize = 9, marker = :xcross)
        end
    end

    Colorbar(fig[1:3, 3]; colormap = colormap, limits = (μlo, μhi),
             label = "μ", width = 14)

    colsize!(fig.layout, 1, Relative(0.55))
    return fig
end

# ----------------------------------------------------------------------------
# 3.  Scalar diagnostics vs μ  (often more readable than the 3-D picture)
# ----------------------------------------------------------------------------

"""
    plot_ellipsoid_summaries(Xs, mus; conf = 0.95)

Centre, principal semi-axis lengths and volume of the ellipsoid as functions
of μ.
"""
function plot_ellipsoid_summaries(Xs, mus; conf = 0.95,
                                  labels = ("θ₁", "θ₂", "θ₃"))
    ells = [fit_ellipsoid(X; conf) for X in Xs]
    fig  = Figure(size = (1050, 320))

    ax1 = Axis(fig[1, 1]; xlabel = "μ", ylabel = "centre")
    for k in 1:3
        y = [e[1][k] for e in ells]
        lines!(ax1, mus, y; label = labels[k]); scatter!(ax1, mus, y)
    end
    axislegend(ax1; position = :lt)

    ax2 = Axis(fig[1, 2]; xlabel = "μ", ylabel = "semi-axis length")
    for k in 1:3
        y = [e[3] * sqrt(sort(eigvals(Symmetric(e[2])); rev = true)[k]) for e in ells]
        lines!(ax2, mus, y; label = "a$k"); scatter!(ax2, mus, y)
    end
    axislegend(ax2; position = :lt)

    ax3 = Axis(fig[1, 3]; xlabel = "μ", ylabel = "volume")
    v = [(4/3)π * e[3]^3 * sqrt(det(Symmetric(e[2]))) for e in ells]
    lines!(ax3, mus, v); scatter!(ax3, mus, v)

    return fig
end

# ----------------------------------------------------------------------------
# 4.  Synthetic data so the script runs out of the box — delete once you plug
#     in your own solutions.
# ----------------------------------------------------------------------------

function demo_data(; mus = range(0, 1; length = 6), n = 100, seed = 42)
    rng = MersenneTwister(seed)
    Xs = Matrix{Float64}[]
    for μ in mus
        c = [1.0 + 0.6μ, -0.8 + 0.5μ, 0.25 - 0.35μ]     # drifting centre
        α = μ * π/3                                      # rotating covariance
        R = [cos(α) -sin(α) 0.0; sin(α) cos(α) 0.0; 0.0 0.0 1.0]
        S = Diagonal([0.30 + 0.10μ, 0.12, 0.05 + 0.06μ]) # anisotropic spread
        L = R * S
        push!(Xs, c' .+ randn(rng, n, 3) * L')
    end
    return Xs, collect(mus)
end
