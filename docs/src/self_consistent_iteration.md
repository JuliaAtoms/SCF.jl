# Self-Consistent Iteration

```@meta
CurrentModule = SCF
```

```@docs
scf!
scf_iteration!
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
mul!
```
