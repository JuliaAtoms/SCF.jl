"""
    solve_orbital_equation!(Pj, eq, method, tol)

Solve the orbital equation `eq` and store the result in the radial
orbital vector of expansion coefficients `Pj`. The solution is
computed using `method`; valid choices are

1. `:arnoldi`, which requires that `hamiltonian(eq)` supports
   [`KrylovWrapper`](@ref), and

2. `:arnoldi_shift_invert`, which iterates `(H-σ*I)⁻¹` and requires
   that `hamiltonian(eq)` supports [`KrylovWrapper`](@ref) *and*
   provides an overload for `MatrixFactorizations.preconditioner`. The
   shift is automatically chosen as `1.1ϵ` where `ϵ` is the (current
   estimate of the) orbital energy of the orbital governed by `eq`.

Both methods are controlled by the stopping tolerance `tol`.
"""
function solve_orbital_equation!(Pj::V, eq::Equation,
                                 method::Symbol, tol;
                                 facttol=√(eps(real(eltype(Pj)))),
                                 io=stdout, verbosity=0,
                                 kwargs...) where {V<:AbstractVector,Equation}
    if method == :arnoldi || method == :arnoldi_shift_invert
        h = hamiltonian(eq)

        schur,history = if method == :arnoldi
            # It would be preferrable if the Arnoldi state could
            # be preserved between iterations, pending
            # https://github.com/haampie/ArnoldiMethod.jl/issues/91
            partialschur(KrylovWrapper(h), nev=1, tol=tol, which=SR())
        elseif method==:arnoldi_shift_invert
            σ = 1.1energy(eq)
            verbosity > 3 && println(io, "Shift = orbital energy of $(eq.orbital) = $(σ) Ha")
            # This is a fantastic waste of memory
            partialschur(ShiftInvert(factorization(KrylovWrapper(h - σ*I);
                                                   tol=facttol,
                                                   isposdefA=true,
                                                   kwargs...), σ),
                         nev=1, tol=tol, which=LR())
        end
        copyto!(Pj, schur.Q[:,1])

        verbosity > 2 && println(io,"Orbital improvement: ", history)
    else
        throw(ArgumentError("Unknown diagonalization method $(method)"))
    end
end
