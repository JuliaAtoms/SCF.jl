"""
    AbstractQuantumSystem

Any implementation of the `AbstractQuantumSystem` interface, used for
self-consistent field calculations, must implement the functions
[`coefficients`](@ref), [`orbitals`](@ref), [`diff`](@ref),
[`normalize!`](@ref).
"""
abstract type AbstractQuantumSystem end

"""
    coefficients(quantum_system)

Retrieves the mixing coefficients of `quantum_system`. Must return a
`view` that the [`scf!`](@ref) routine can modify. To be overloaded by
the implementation of [`AbstractQuantumSystem`](@ref).
"""
coefficients(::QS) where {QS<:AbstractQuantumSystem} =
    throw(ArgumentError("`coefficients` not implemented for $(QS)"))

"""
    orbitals(quantum_system)

Retrieves the orbitals of `quantum_system`. Must return a
`view` that the [`scf!`](@ref) routine can modify. To be overloaded by
the implementation of [`AbstractQuantumSystem`](@ref).
"""
orbitals(::QS) where {QS<:AbstractQuantumSystem} =
    throw(ArgumentError("`orbitals` not implemented for $(QS)"))

"""
    diff(quantum_system[; kwargs...])

Varies the the `quantum_system` with respect to all orbitals. Used to
derive the multi-configurational Hartree–Fock equations. To be overloaded by
the implementation of [`AbstractQuantumSystem`](@ref).
"""
Base.diff(::QS; kwargs...) where {QS<:AbstractQuantumSystem} =
    throw(ArgumentError("`diff` not implemented for $(QS)"))

"""
    normalize!(quantum_system, v)

Normalize the orbital `v` of `quantum_system`. To be overloaded by the
implementation of [`AbstractQuantumSystem`](@ref).
"""
LinearAlgebra.normalize!(::QS, ::AbstractVector) where {QS<:AbstractQuantumSystem} =
    throw(ArgumentError("`normalize!` not implemented for `$QS`"))

export AbstractQuantumSystem
