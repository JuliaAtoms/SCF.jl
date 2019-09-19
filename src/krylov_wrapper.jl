"""
    KrylovWrapper(hamiltonian)

Proxy object used in the Krylov iterations, during orbital
improvement. This is useful, since `hamiltonian` may be defined to act
on objects such as vectors living in function spaces (as
e.g. implemented using ContinuumArrays.jl), whereas the SCF iterations
act on the coefficients directly.
"""
mutable struct KrylovWrapper{T,Hamiltonian}
    hamiltonian::Hamiltonian
end

Base.eltype(A::KrylovWrapper{T}) where T = T
Base.size(A::KrylovWrapper{T,<:AbstractMatrix}) where T = size(A.hamiltonian)
"""
    size(::KrylovWrapper)

Returns the dimension of the `KrylovWrapper`. For Hamiltonians which
are not `<:AbstractMatrix`, this needs to be overloaded.
"""
Base.size(A::KrylovWrapper) = size(A.hamiltonian, A)
Base.size(A::KrylovWrapper, i) = size(A)[i]

IterativeFactorizations.factorization(A::KrylovWrapper; kwargs...) =
    IterativeFactorization(A; kwargs...)

IterativeFactorizations.preconditioner(A::KrylovWrapper) =
    IterativeFactorizations.preconditioner(A.hamiltonian)

"""
    mul!(y, ::KrylovWrapper, x)

Compute the action of the wrapped Hamiltonian on `x` and store it in
`y`. For Hamiltonians which are not `<:AbstractMatrix`, this needs to
be overloaded.
"""
LinearAlgebra.mul!(y::V₁, A::KrylovWrapper{T,Hamiltonian}, x::V₂) where {V₁,V₂,T,Hamiltonian<:AbstractMatrix} =
    mul!(y, A.hamiltonian, x)

Base.show(io::IO, kw::KrylovWrapper{T}) where T =
    write(io, "KrylovWrapper{$T} of size $(size(kw))")

"""
    OrthogonalKrylovWrapper(h,A,x′,C,S,tmp)

Proxy object to perform Krylov iterations of the Hamiltonian `h`,
represented by the [`KrylovWrapper`](@ref) `A`, in the subspace
orthogonal to `C`. The component of the test vector orthogonal to `C`
is stored in `x′`, before acting with `A`, and finally orthogonalized
again, i.e. `y = (1-C*C')*A*(1-C*C')*x = (1-C*C')*A*x′`. `S` is the
overlap matrix of the basis functions, used to compute the inner
products.
"""
struct OrthogonalKrylovWrapper{H,K<:KrylovWrapper,V,M,SM}
    h::H
    A::K
    x′::V
    C::M
    S::SM
    SC::M
end

function OrthogonalKrylovWrapper(h, A, nC, S)
    m = size(A, 1)
    T = eltype(A)
    x′ = zeros(T, m)
    C = zeros(T, m, nC)
    SC = zeros(T, m, nC)
    OrthogonalKrylovWrapper(h, A, x′, C, S, SC)
end

Base.eltype(A::OrthogonalKrylovWrapper) = eltype(A.A)
Base.size(A::OrthogonalKrylovWrapper) = size(A.A)
Base.size(A::OrthogonalKrylovWrapper, i) = size(A.A, i)

function orthogonalize!(x′, C, SC, x)
    for j = 1:size(C,2)
        c = dot(view(SC, :, j), x)
        BLAS.axpy!(-c, view(C, :, j), x′)
    end
end

function LinearAlgebra.mul!(y, A::OrthogonalKrylovWrapper, x)
    # We cannot simply project `x` into the subspace orthogonal to
    # `C`, since `x` is a test vector from e.g. a Krylov iteration,
    # i.e. it is part of the Krylov space and must not be
    # modified. Hence, we make our own copy, and project out the
    # components parallel to `C` from that one instead.
    copyto!(A.x′, x)
    # This step is necessary for non-orthogonal bases, since the inner
    # product is defined as `a'S*b`, where `S=B'B`, which is
    # non-diagonal for non-orthogonal bases.
    mul!(A.SC, A.S, A.C)

    # x′ <- (1 - C*C')*x
    orthogonalize!(A.x′, A.C, A.SC, x)
    # y <- A*x′
    mul!(reshape(y, :), A.A, A.x′)
    # y <- (1 - C*C')*y
    orthogonalize!(y, A.C, A.SC, y)
end

Base.show(io::IO, okw::OrthogonalKrylovWrapper) =
    write(io, "OrthogonalKrylovWrapper of size $(size(okw)), orthogonal against subspace of $(size(okw.C,2)) vectors")
