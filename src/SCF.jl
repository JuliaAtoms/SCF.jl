module SCF

using LinearAlgebra

using SolverTraces
import SolverTraces: base_exp

using UnicodeFun
using Formatting

if VERSION < v"1.1-DEV"
    isnothing(::Nothing) = true
    isnothing(::Any) = false
end

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
    show(io, "text/plain", fock.equations)
end

Base.view(::Fock{Q}, args...) where Q =
    throw(ArgumentError("`view` not implemented for `Fock{$Q}`"))

LinearAlgebra.ldiv!(::Fock{Q}, v) where Q =
    throw(ArgumentError("`ldiv!` not implemented for `Fock{$Q}`"))

"""
    scf!([fun!, ]fock[; ω=0, max_iter=200, tol=1e-8, verbosity=1])

Perform the SCF iterations for the `fock` operator. The actual
solution of the secular equations is performed by the implementation
of `ldiv!(y, fock, x)`. The optional `fun!` argument allows for extra
operations to be performed every SCF cycle (such as rotations, etc).

`ω` is a relaxation parameter; the coefficients are updated as `cᵢ₊₁ =
(1-ω)c̃ + ωcᵢ` where `c̃` is the solution to secular equations in the
current iteration and `cᵢ` the previous set of coefficients. The
default (`ω=0`) is to only use this solution.

The SCF procedure continues until either the amount of iterations
equals `max_iter` or the change in the coefficients is below `tol`.
"""
function scf!(fun!::Function, fock::Fock{Q};
              ω=0, max_iter=200, tol=1e-8,
              verbosity=1) where Q
    trace,tolerance = if verbosity > 1
        trace = SolverTrace(max_iter,
                            CurrentStep(max_iter,
                                        lc=LinearColorant(max_iter,1,SolverTraces.red_green_scale()),
                                        header="Iteration"))

        tolerance = Tolerance(tol, print_target=false)
        push!(trace, tolerance)
        trace,tolerance
    else
        nothing,nothing
    end

    orbitals = view(fock, :, :)
    norb = size(orbitals,2)
    coeffs = copy(orbitals)
    c̃ = similar(coeffs)

    Δ = [Inf for j in 1:norb]

    if verbosity ≥ 2
        println("Self-consistent field calculation of")
        print("- ")
        display(fock.quantum_system)
        println("- Maximum amount of iterations: $(max_iter)")
        tb,te = base_exp(tol)
        println("- Stopping tolerance: ", format(tolerance.tol_fmt, tb, to_superscript(te)))
        ω != 0 && println("- Successive relaxation: ω = $(ω)")
        println()
    end

    print_header(trace)
    for i = 1:max_iter
        ldiv!(fock, coeffs)
        fun!(coeffs)

        for j = 1:norb
            Δ[j] = 1.0 - dot(coeffs[:,j],orbitals[:,j])
        end

        aΔ = norm(Δ)
        if aΔ < tol
            println()
            break
        end
        isnothing(tolerance) || (tolerance.current = aΔ)

        if ω == 0
            copyto!(orbitals, coeffs)
        else
            # Relaxation; see e.g. p. 490 of
            #
            # - Fischer, C. F., & Guo, W. (1990). Spline algorithms for the
            #   hartree-fock equation for the helium ground state. Journal of
            #   Computational Physics, 90(2),
            #   486–496. http://dx.doi.org/10.1016/0021-9991(90)90176-2
            lmul!(ω, orbitals)
            lmul!(1-ω, coeffs)
            orbitals[:] += coeffs[:] # Is this efficient?
        end

        SolverTraces.next!(trace)
    end

    norm(Δ) > tol && @warn "Desired tolerance $(tol) not reached in $(max_iter) iterations"

    fock
end

scf!(fock::Fock{Q}; kwargs...) where Q =
    scf!((_)->nothing, fock; kwargs...)

export Fock, scf!

end # module
