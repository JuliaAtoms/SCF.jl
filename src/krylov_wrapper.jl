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

MatrixFactorizations.factorization(A::KrylovWrapper; kwargs...) =
    IterativeFactorization(A; kwargs...)

MatrixFactorizations.preconditioner(A::KrylovWrapper) =
    MatrixFactorizations.preconditioner(A.hamiltonian)

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
