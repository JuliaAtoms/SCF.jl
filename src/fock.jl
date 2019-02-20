"""
    Fock(quantum_system, equations)

A Fock operator consists of a `quantum_system`, from which `equations`
are variationally derived. `equations` must be iterable, where each
element corresponds to the equation for one orbital, and must provide
overloads for the functions [`energy`](@ref),
[`energy_matrix!`](@ref), and [`hamiltonian`](@ref).  Additionally,
[`update!`](@ref) must be provided for `equations`, to prepare the
equation system for the next iteration.
"""
mutable struct Fock{Q<:AbstractQuantumSystem,Equations}
    quantum_system::Q
    equations::Equations
end

Fock(quantum_system::Q; kwargs...) where Q =
    Fock(quantum_system, diff(quantum_system; kwargs...))

function Base.show(io::IO, ::MIME"text/plain", fock::Fock{Q,E}) where {Q,E}
    write(io, "Fock operator with\n- quantum system: ")
    show(io, "text/plain", fock.quantum_system)
    write(io, "\n- SCF equations:  ")
    for eq in fock.equations
        write(io, "\n  - ")
        show(io, "text/plain", eq)
    end
end

Base.view(::Fock{Q}, args...) where Q =
    throw(ArgumentError("`view` not implemented for `Fock{$Q}`"))

"""
    energy(equation::Equation)

Computes the orbital energy of `equation`. _To be overloaded by the
user._
"""
energy(::Equation) where Equation =
    throw(ArgumentError("`energy` not implemented for `$Equation`"))
"""
    energy(fock::Fock)

Calculates the total energy of the system by summing the orbital energies.
"""
energy(fock::Fock) = sum(energy(eq) for eq in fock.equations)

"""
    energy_matrix!(H::AbstractMatrix, equation::Equation)

Assemble the energy matrix for `equation`, i.e. the matrix where the
elements are the energies corresponding to the mixing of
configurations, basically the Hamiltonian matrix of the corresponding
orbital. This matrix, sandwiched between the mixing coefficients,
gives the orbital energy. _To be overloaded by the user._ NB the
matrix elements should be added to, i.e. *not overwritten*, rather,
the overload of `energy_matrix!` should compute (the equivalent of) `H
+= H_i`, where `H_i` is the energy matrix for `equation`.
"""
energy_matrix!(H::AbstractMatrix,::Equation) where Equation =
    throw(ArgumentError("`energy_matrix!` not implemented for `$Equation`"))

"""
    energy_matrix!(H::AbstractMatrix, fock::Fock)

Calculates the total energy matrix of the system by summing the energy
matrices of the different orbital equations of `fock`. This overwrites
the entries of `H`.
"""
function energy_matrix!(H::HM,fock::F) where {HM<:AbstractMatrix,F<:Fock}
    H .= zero(eltype(H))
    foreach(eq -> energy_matrix!(H, eq), fock.equations)
    H
end

"""
    hamiltonian(equation::Equation)

Returns the orbital Hamiltonian of `equation`. _To be overloaded by the
user._
"""
hamiltonian(equation::Equation) where Equation =
    throw(ArgumentError("`hamiltonian` not implemented for `$Equation`"))

"""
    update!(eqs; kwargs...)

Update the equation system `eqs` for the current iteration. _To be
overloaded by the user._
"""
update!(eqs; kwargs...) =
    throw(ArgumentError("`update!` not implemented for $(typeof(eqs))"))

"""
    rotate_max_lobe!(v)

Rotate the vector `v` such that the largest lobe has positive sign.
"""
function rotate_max_lobe!(v::V) where {V<:AbstractVector}
    i = argmax(abs.(v))
    lmul!(sign(v[i]), v)
    v
end

"""
    norm_rot!(fock, v)

Normalize and rotate the eigenvector `v` such that the largest
lobe has positive sign.
"""
function norm_rot!(fock::Fock, v::V) where {V<:AbstractVector}
    normalize!(fock.quantum_system, v)
    rotate_max_lobe!(v)
end
