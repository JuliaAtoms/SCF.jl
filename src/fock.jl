"""
    Fock(quantum_system, equations)

A Fock operator consists of a `quantum_system`, from which `equations`
are variationally derived.
"""
mutable struct Fock{Q<:AbstractQuantumSystem,E}
    quantum_system::Q
    equations::Vector{E}
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

energy(fock::Fock) = sum(energy(eq) for eq in fock.equations)

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

update!(eq; kwargs...) =
    throw(ArgumentError("`update!` not implemented for $(typeof(eq))"))
