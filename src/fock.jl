mutable struct Fock{Q,E}
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

LinearAlgebra.normalize!(::Fock{Q}, ::AbstractVector) where Q =
    throw(ArgumentError("`normalize!` not implemented for `Fock{$Q}`"))

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
    normalize!(fock, v)
    rotate_max_lobe!(v)
end

update!(eq; kwargs...) =
    throw(ArgumentError("`update!` not implemented for $(typeof(eq))"))

function LinearAlgebra.ldiv!(fock::Fock{Q,E}, c::M;
                             verbosity=0, method=:arnoldi,
                             tol=1e-10, kwargs...) where {Q,E,M<:AbstractVecOrMat}
    verbosity > 0 && println("Solving secular equations using $(method)")
    for (j,eq) in enumerate(fock.equations)
        print_block() do io
            # Ideally, all direct potentials should be shared among the
            # equations, such that they are only calculated once per
            # iteration.
            verbosity > 1 && println(io, eq)
            update!(eq; verbosity=max(0,verbosity-2), io=io)

            vcj = view(c,:,j)

            if method==:arnoldi
                # It would be preferrable if the Arnoldi state could
                # be preserved between iterations, pending
                # https://github.com/haampie/ArnoldiMethod.jl/issues/91
                K = KrylovWrapper(eq)
                schur,history = partialschur(K, nev=1, tol=tol, which=SR())
                copyto!(vcj, schur.Q[:,1])
                #= elseif method==:arnoldi_shift_invert
                # If we figure out a method of factorizing the HF
                # Hamiltonian at not a too large cost, we could use
                # the shift-and-invert trick here as well, using the
                # "hydrogenic" (i.e. without repulsion integrals)
                # energy as a lower bound.
                =#
            else
                throw(ArgumentError("Unknown diagonalization method $(method)"))
            end

            # Normalize eigenvector and rotate it such that the
            # largest lobe is positive.
            norm_rot!(fock, vcj)

            if verbosity > 2
                println(io, "Secular equation solution: ", history)
                println(io, "Change in equation $j: ",
                        norm(c[:,j] - eq.ϕ.mul.factors[2]))
            end
        end
    end
    fock
end
