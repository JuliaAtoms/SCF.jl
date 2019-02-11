"""
    solve_secular_problem!(H, c, fock)

Form the energy matrix, store it in `H`, and then solve the secular
problem `Hc = Ec` for the lowest eigenvalue.
"""
function solve_secular_problem!(H::M, c::C, fock::F;
                                tol=√(eps()),
                                verbosity=0) where {T,M<:AbstractMatrix{T},
                                                    C<:AbstractVector{T}, F<:Fock}
    verbosity > 0 && println("Solving secular problem")

    energy_matrix!(H, fock)

    # This could be more efficient if we could use c as the initial
    # guess for the Arnoldi procedure.
    schur,history = partialschur(H, nev=1, tol=tol, which=SR())
    copyto!(c, schur.Q[:,1])

    normalize!(c)
end
