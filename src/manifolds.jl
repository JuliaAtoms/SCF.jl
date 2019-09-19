import Optim: retract!, project_tangent!

# * Sphere
metric_dot(a, b, ::UniformScaling{Bool}) = dot(a,b)
metric_dot(a, b, g::UniformScaling) = dot(a,b)*g.λ
metric_dot(a, b, g) = dot(a,g*b)

metric_norm(x, g::UniformScaling) = norm(x)*√(g.λ)
metric_norm(x, g) = √(metric_dot(x,x,g))

metric_normalize!(x, ::UniformScaling{Bool}) = normalize!(x)
function metric_normalize!(x, g)
    x ./= metric_norm(x,g)
end

"""Spherical manifold {|x| = 1}."""
struct MetricSphere{M} <: Manifold
    "Metric used for inner products, i.e. `⟨a|b⟩_g = gₘₙaᵐbⁿ"
    g::M
end
MetricSphere() = MetricSphere(I)
Optim.retract!(S::MetricSphere, x) = metric_normalize!(x, S.g)
Optim.project_tangent!(S::MetricSphere,g,x) = (g .-= real(metric_dot(x, g, S.g)).*x)

# * Stiefel

matfun(f::Function, ee) = ee.vectors*Diagonal(f.(ee.values))*ee.vectors'

# For orthogonal bases, e.g. FEDVRQuasi.jl
löwdin_transform(g::UniformScaling{Bool}) = g,g
# For quasi-orthgonal bases, e.g. FiniteDifferencesQuasi.jl
löwdin_transform(g::UniformScaling) = inv(√(g.λ))*I,√(g.λ)*I
# For non-orthogonal bases, e.g. BSplinesQuasi.jl
function löwdin_transform(g::AbstractMatrix)
    ee = eigen(g)
    matfun(λ -> inv(√(λ)), ee), matfun(λ -> √(λ), ee)
end

löwdin_transform!(v, ::UniformScaling{Bool}) = v
löwdin_transform!(v, S⁻ᴴ) = (v .= S⁻ᴴ*v)

abstract type MetricStiefel <: Stiefel end
struct MetricStiefel_CholQR{L} <: MetricStiefel
    "Löwdin transform to orthogonal basis"
    S⁻ᴴ::L
    "Löwdin transform to non-orthogonal basis"
    Sᴴ::L
end
struct MetricStiefel_SVD{L} <: MetricStiefel
    "Löwdin transform to orthogonal basis"
    S⁻ᴴ::L
    "Löwdin transform to non-orthogonal basis"
    Sᴴ::L
end
function MetricStiefel(g, retraction=:SVD)
    S⁻ᴴ,Sᴴ = löwdin_transform(g)
    if retraction == :CholQR
        MetricStiefel_CholQR(S⁻ᴴ,Sᴴ)
    elseif retraction == :SVD
        MetricStiefel_SVD(S⁻ᴴ,Sᴴ)
    else
        throw(ArgumentError("Unknown retraction $retraction"))
    end
end
function retract!(St::MetricStiefel_SVD, X)
    löwdin_transform!(X, St.Sᴴ)
    U,S,V = svd(X)
    X .= U*V'
    löwdin_transform!(X, St.S⁻ᴴ)
end
function retract!(St::MetricStiefel_CholQR, X)
    löwdin_transform!(X, St.Sᴴ)
    overlap = X'X
    X .= X/cholesky(overlap).L
    löwdin_transform!(X, St.S⁻ᴴ)
end
# For functions depending only on the subspace spanned by X, we always
# have G = A*X for some A, and so X'G = G'X, and Stiefel == Grassmann
# Edelman et al. have G .-= X*G'X (2.53), corresponding to a different
# metric ("canonical metric"). We follow Absil et al. here and use the
# metric inherited from Nxn matrices.
function project_tangent!(St::MetricStiefel, G, X)
    löwdin_transform!(X, St.Sᴴ)
    löwdin_transform!(G, St.Sᴴ)
    XG = X'G
    G .-= X*((XG .+ XG')./2)
    löwdin_transform!(X, St.S⁻ᴴ)
    löwdin_transform!(G, St.S⁻ᴴ)
end

# * Manifold setup

function setup_manifold(fock::Fock, m, n)
    # all(isone, length.(fock.symmetries)) || 
    #     @warn "This is only valid for system without any orthogonality constraints"
    # Optim.PowerManifold(MetricSphere(fock.S), (m,), (n,))
    MetricStiefel(fock.S)
end
