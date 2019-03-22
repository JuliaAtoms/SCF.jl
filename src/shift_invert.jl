"""
    ShiftInvert(A⁻¹, σ)

Structure representing the inverse of a shifted matrix `Ã = A +
σ*I`. This is used for finding eigenpairs of `Ã` in the neighbourhood
of `σ`; useful in power methods such as the Krylov iteration, which
converge the convex hull of the eigenspectrum first.
"""
struct ShiftInvert{Factorization, T}
    A⁻¹::Factorization
    σ::T
end

Base.size(S::ShiftInvert, args...) = size(S.A⁻¹, args...)
Base.eltype(S::ShiftInvert) = eltype(S.A⁻¹)

# """
#     ShiftInvert(A⁻¹, σ)

# Construct the [`ShiftInvert`](@ref) structure from an already computed
# factorization `A⁻¹` of the original matrix shifted by `σ`.
# """
# ShiftInvert(A⁻¹::Factorization, σ::T) where {Factorization, T} =
#     ShiftInvert{Factorization,T}(A⁻¹, σ)

"""
    ShiftInvert(A[, σ])

Construct the [`ShiftInvert`](@ref) structure from the factorization
of the matrix `A` shifted by `σ`. If no shift is specified, the lower
bound as calculated by [`gershgorin_bounds`](@ref) is used by default.
"""
ShiftInvert(A::M, σ::T=leftendpoint(gershgorin_bounds(A));
            factorization=factorize) where {M<:AbstractMatrix, T} =
    ShiftInvert(factorization(A - σ*I), σ)

"""
    mul!(y, si::ShiftInvert, x)

Compute the action of the shifted-and-inverted matrix `si.A⁻¹` on `x`
and store the result in `y`. It is assumed that the factorization does
not depend on the initial contents of `y`; if that is not the case, an
overload has to be provided.
"""
LinearAlgebra.mul!(y, si::ShiftInvert, x) =
    ldiv!(y, si.A⁻¹, x)

"""
    mul!(y, si::ShiftInvert{<:IterativeFactorization}, x)

Compute the action of the shifted-and-inverted matrix `si.A⁻¹` on `x`
and store the result in `y` for factorizations based on
[IterativeSolvers.jl](https://github.com/JuliaMath/IterativeSolvers.jl).
"""
function LinearAlgebra.mul!(y, si::ShiftInvert{<:IterativeFactorization}, x)
    y .= false # Strong zero to get rid of NaNs
    ldiv!(y, si.A⁻¹, x)
end

export ShiftInvert
