# Fock Operators

```@meta
CurrentModule = SCF
```

The Fock operator consists of an [`AbstractQuantumSystem`](@ref) and a
set of coupled integro-differential equations, the solution of which
is the objective of the self-consistent field procedure. For the
solution process to work, the set of equations need to implement a few
methods:

1. It must fulfil Julia's iteration interface, i.e. each element must
   be the equation for a single orbital, which is solved independently
   from the other equations, but with the other orbitals as inputs.
2. [`energy_matrix!`](@ref) which calculates the energy matrix
   $\mat{H}_i$ for orbital equation $i$, where
   $\vec{c}^H\mat{H}_i\vec{c}$ gives the orbital energy for the
   correspond orbital. The overall energy matrix
   $\mat{H}=\sum_i\mat{H}_i$ is used to solve the secular problem for
   the mixing coefficients.
3. [`hamiltonian`](@ref) which returns the Hamiltonian corresponding
   to one orbital equation.
4. [`update!`](@ref) which recomputes all orbital-dependent integrals,
   shared among the equations of the equation system.

```@docs
Fock
norm_rot!
rotate_first_lobe!
energy_matrix!
hamiltonian
update!
```

