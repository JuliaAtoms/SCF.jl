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

## Shift-and-invert

To improve the rate of convergence of the Krylov iterations, instead
of computing the eigenvectors of $\Hamiltonian$ when improving the
orbitals, we can use the shift-and-invert strategy and instead iterate
$(\Hamiltonian-\sigma I)^{-1}$ with a shift $\sigma$ that is located
slightly below the target eigenvalue. The eigenvectors are shared
between $\Hamiltonian$ and $(\Hamiltonian-\sigma I)^{-1}$, but
convergence towards the low-lying states is faster for
$(\Hamiltonian-\sigma I)^{-1}$.

This approach is complicated by the fact that the orbital Hamiltonian
$\Hamiltonian$ is an integro-differential operator, which furthermore
depends on the other orbitals through direct, exchange, and
configuration interaction, and thus $\Hamiltonian$ is not easily
factorizable. However, we can compute the action of
$(\Hamiltonian-\sigma I)^{-1}$ using an [iterative
solver](https://github.com/JuliaMath/IterativeSolvers.jl) with a
preconditioner constructed from all terms of $\Hamiltonian$ except the
integral operators (exchange interaction) and the source terms
(configuration interaction). As mentioned in the documentation for
[`solve_orbital_equation!`](@ref), the orbital Hamiltonian must
support [`KrylovWrapper`](@ref).

```@docs
ShiftInvert
ShiftInvert(A::M, σ::T=leftendpoint(gershgorin_bounds(A)); factorization=factorize) where {M<:AbstractMatrix, T}
LinearAlgebra.mul!(y, si::ShiftInvert, x)
```

### Gershgorin's circle theorem
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
