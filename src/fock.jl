"""
    Fock(quantum_system, equations)

A Fock operator consists of a `quantum_system`, from which `equations`
are variationally derived. `equations` must provide an overload for
[`energy_matrix!`](@ref). `equations` must be iterable, where each
element corresponds to the equation for one orbital, and must provide
an overload for [`hamiltonian`](@ref) and
[`energy`](@ref). Additionally, [`update!`](@ref) must be provided for
`equations`, to prepare the equation system for the next iteration.
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
    energy_matrix!(H::AbstractMatrix, equations[, which=:total_energy])

Calculates the energy matrix of the system of `equations`. This
overwrites the entries of `H`. _To be overloaded by the user; must
support `which=:total_energy`, `which=:double_counted_energy`, and
`which=:kinetic_energy`._
"""
energy_matrix!(H::HM, equations::E, which::Symbol=:total_energy) where {HM<:AbstractMatrix,E} =
    throw(ArgumentError("`energy_matrix!` not implemented for $E"))

"""
    energy_matrix!(H::AbstractMatrix, fock::Fock[, which=:total_energy])

Calculates the energy matrix of the quantum system of `fock`. This
overwrites the entries of `H`.
"""
function energy_matrix!(H::HM, fock::F, which::Symbol=:total_energy) where {HM<:AbstractMatrix,F<:Fock}
    H .= zero(eltype(H))
    energy_matrix!(H, fock.equations, which)
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
    energy(equation::Equation[, term=:total_energy])

Returns the orbital energy for the orbital governed by `equation`. _To
be overloaded by the user._
"""
energy(equation::Equation, which::Symbol=:total_energy) where Equation =
    throw(ArgumentError("`energy` not implemented for `$Equation`"))

"""
    update!(eqs; kwargs...)

Update the equation system `eqs`, for the current iteration. _To be
overloaded by the user._
"""
update!(eqs; kwargs...) =
    throw(ArgumentError("`update!` not implemented for $(typeof(eqs))"))

"""
    update!(eqs, quantum_system; kwargs...)

Update the equation system `eqs` with respect to `quantum_system`, for
the current iteration. _To be overloaded by the user._
"""
update!(eqs, quantum_system; kwargs...) =
    throw(ArgumentError("`update!` not implemented for $(typeof(eqs)), $(typeof(quantum_system))"))

update!(fock::Fock; kwargs...) = update!(fock.equations)

"""
    rotate_first_lobe!(v)

Rotate the vector `v` such that the first lobe has positive sign.
"""
function rotate_first_lobe!(v::V) where {V<:AbstractVector}
    lmul!(sign(v[1]), v)
    v
end

"""
    norm_rot!(fock, v)

Normalize and rotate the eigenvector `v` such that the first
lobe has positive sign.
"""
function norm_rot!(fock::Fock, v::V) where {V<:AbstractVector}
    normalize!(fock.quantum_system, v)
    rotate_first_lobe!(v)
end

"""
    Matrix(fock)

Compute the Fock matrix, where each element corresponds to a matrix
element between two orbitals.
"""
function LinearAlgebra.Matrix(fock::Fock)
    m = length(fock.equations)
    F = Matrix{Float64}(undef, m, m)

    orbs = orbitals(fock.quantum_system)
    tmp = Vector{eltype(orbs)}(undef, size(orbs, 1))

    for (i,eqi) in enumerate(fock.equations)
        hi = KrylovWrapper(hamiltonian(eqi))
        orbi = view(orbs, :, i)

        for (j,eqj) in enumerate(fock.equations)
            mul!(tmp, hi, orbi)
            F[i,j] = real(dot(view(orbs, :, j), tmp))
        end
    end
    F
end
