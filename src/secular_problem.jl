"""
    solve_secular_problem!(H, c, fock)

Form the energy matrix, store it in `H`, and then solve the secular
problem `Hc = Ec` for the lowest eigenvalue.
"""
function solve_secular_problem!(H::M, c::C, fock::F, update_mixing_coefficients;
                                tol=√(eps()),
                                verbosity=0) where {T,M<:AbstractMatrix{T},
                                                    C<:AbstractVector{T}, F<:Fock}
    energy_matrix!(H, fock)

    if update_mixing_coefficients
        if length(c) == 1
            c[1] = one(T)
            return c
        end

        verbosity > 0 && println("Solving secular problem")

        # This could be more efficient if we could use c as the initial
        # guess for the Arnoldi procedure.
        schur,history = partialschur(H, nev=1, tol=tol, which=SR())
        copyto!(c, schur.Q[:,1])
    end

    normalize!(c)
end
