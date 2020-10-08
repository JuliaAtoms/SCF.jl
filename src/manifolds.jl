import Optim: retract!, project_tangent!

# * Sphere

metric_norm(x, g::UniformScaling) = norm(x)*√(g.λ)
metric_norm(x, g) = √(dot(x,g,x))

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
Optim.project_tangent!(S::MetricSphere,g,x) = (g .-= real(dot(x, S.g, g)).*x)

# * Stiefel

matfun(f::Function, ee) = ee.vectors*Diagonal(f.(ee.values))*ee.vectors'

# For orthogonal bases, e.g. FE-DVR
löwdin_transform(g::UniformScaling{Bool}) = g,g
# For quasi-orthogonal bases, e.g. finite-differences
löwdin_transform(g::UniformScaling) = inv(√(g.λ))*I,√(g.λ)*I
löwdin_transform(g::Diagonal) = inv(√(g)),√(g)
# For non-orthogonal bases, e.g. B-splines
function löwdin_transform(g::AbstractMatrix)
    ee = eigen(g)
    matfun(λ -> inv(√(λ)), ee), matfun(λ -> √(λ), ee)
end
function löwdin_transform(g::BandedMatrix)
    m,n = size(g)
    m == n || throw(DimensionMismatch("Löwdin transform only valid for square matrices"))
    l,u = bandwidths(g)
    if l == u == 0
        löwdin_transform(Diagonal(g))
    else
        löwdin_transform(Symmetric(g))
    end
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

# * Sphere & Stiefel combo

struct SphereStiefelCombo{Sp,St,M,F} <: Manifold
    "Indices of columns that only need normalization"
    normal_indidices::Vector{Int}
    "Manifold used to normalize columns"
    sphere::Sp
    "Sets of indices of columns that should be mutually orthogonalized"
    ortho_indidices::Vector{Vector{Int}}
    "Manifold used to orthogonalize sets of columns"
    stiefel::St
    "Temporary storage for orthgonalization"
    ortho_X::M
    "Temporary storage for orthgonalization"
    ortho_G::M
    f::F
end

function copytotmp!(dst, src, ortho_indices)
    for (i,j) in enumerate(ortho_indices)
        copyto!(view(dst, :, i), view(src, :, j))
    end
end

function copyfromtmp!(dst, src, ortho_indices)
    for (i,j) in enumerate(ortho_indices)
        copyto!(view(dst, :, j), view(src, :, i))
    end
end

function retract!(ssc::SphereStiefelCombo, X)
    retract!(ssc.sphere, view(X, :, ssc.normal_indidices))
    for o in ssc.ortho_indidices
        copytotmp!(ssc.ortho_X, X, o)
        retract!(ssc.stiefel, view(ssc.ortho_X, :, 1:length(o)))
        copyfromtmp!(X, ssc.ortho_X, o)
    end
    X
end

function project_tangent!(ssc::SphereStiefelCombo, G, X)
    project_tangent!(ssc.sphere, view(G, :, ssc.normal_indidices), view(X, :, ssc.normal_indidices))
    for o in ssc.ortho_indidices
        sel = 1:length(o)
        copytotmp!(ssc.ortho_G, G, o)
        copytotmp!(ssc.ortho_X, X, o)
        project_tangent!(ssc.stiefel, view(ssc.ortho_G, :, sel), view(ssc.ortho_X, :, sel))
        copyfromtmp!(G, ssc.ortho_G, o)
        copyfromtmp!(X, ssc.ortho_X, o)
    end
    G
end

# * Manifold setup

function setup_manifold(f::FockProblem)
    fock = f.fock
    P = orbitals(fock.quantum_system)
    m,n = size(P)

    normal_indices = sort(reduce(vcat, filter(s -> length(s)==1, fock.symmetries), init=Int[]))
    ortho_indices = filter(s -> length(s)≠1, fock.symmetries)
    n_ortho = length.(ortho_indices)

    length(normal_indices) + sum(n_ortho) == n ||
        throw(DimensionMismatch("Number of orbitals does not agree with amount of symmetry specified"))

    sphere = Optim.PowerManifold(MetricSphere(fock.S), (m,), (length(normal_indices),))
    stiefel = MetricStiefel(fock.S)

    nn = isempty(n_ortho) ? 0 : maximum(n_ortho)

    T = eltype(P)
    ortho_X = Matrix{T}(undef, m, nn)
    ortho_G = Matrix{T}(undef, m, nn)

    SphereStiefelCombo(normal_indices, sphere,
                       ortho_indices, stiefel,
                       ortho_X, ortho_G, f)
end
