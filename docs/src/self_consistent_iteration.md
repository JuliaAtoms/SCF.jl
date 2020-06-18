# Self-Consistent Iteration

```@meta
CurrentModule = SCF
```

```@docs
scf!
scf_iteration!
solve_orbital_equation!
solve_secular_problem!
```

## KrylovWrapper

The orbital improvements are performed via diagonalization of the
integro-differential equations for each orbital. The diagonalization
procedure is the Arnoldi method, which requires the repeated action of
the Hamiltonian on a vector. The [`KrylovWrapper`](@ref) object is
used to wrap Hamiltonians which are not simple `<:AbstractMatrix`
objects, but which require a specialized implementation of
[`mul!`](@ref) to act on a `<:AbstractVector` of coefficients. An
example is the kind of Hamiltonians implemented in the
[Atoms.jl](https://github.com/JuliaAtoms/Atoms.jl) library, which are
based on the function space algebra from
[ContinuumArrays.jl](https://github.com/JuliaApproximation/ContinuumArrays.jl).

```@docs
KrylovWrapper
Base.size(::KrylovWrapper)
LinearAlgebra.mul!(y::V₁, A::KrylovWrapper{T,Hamiltonian}, x::V₂) where {V₁,V₂,T,Hamiltonian<:AbstractMatrix}
```

## Gershgorin's circle theorem
[Gershgorin's circle
theorem](https://en.wikipedia.org/wiki/Gershgorin_circle_theorem) can
be used to estimate the bounds on the eigenspectrum of a matrix. It is
always correct, but only a useful estimate when the matrix is
diagonally dominant. Furthermore, [`gershgorin_bounds`](@ref) is
limited to the case of Hermitian matrices, whose eigenvalues lie on
the real line. If an experimental energy is known, that is almost
always a better lower bound to use for the [Shift-and-invert](@ref)
strategy.

```@docs
gershgorin_disc
gershgorin_discs
gershgorin_bounds
```
