"""
    solve_orbital_equation!(Pj, okw, method, tol)

Solve the orbital equation for the `j`th orbital, by solving the
eigenproblem of the [`OrthogonalKrylovWrapper`](@ref) operator `okw`,
and store the result in the radial orbital vector of expansion
coefficients `Pj`. The solution is computed using `method`; valid
choices are

1. `:arnoldi`, which uses ArnoldiMethod.jl. Really small grid
   spacings/large grid extents can be problematic with this method,
   due to high condition numbers.

2. `:lobpcg`, which uses IterativeSolvers.jl. Too coarse grids can
   cause never-ending loops or errors from the Cholesky decomposition,
   complaining about non-positive definite subproblems.

Both methods are controlled by the stopping tolerance `tol`.
"""
function solve_orbital_equation!(Pj, okw, method, tol;
                                 facttol=√(eps(real(eltype(Pj)))),
                                 io=stdout, verbosity=0,
                                 kwargs...)
    verbosity > 0 && println(io, "Improving orbital using $method")

    update!(okw)

    ϕ = if method == :arnoldi
        # It would be preferrable if the Arnoldi state could
        # be preserved between iterations, pending
        # https://github.com/haampie/ArnoldiMethod.jl/issues/91
        schur,history = partialschur(okw, nev=1, tol=tol, which=SR())
        verbosity > 2 && println(io,"Orbital improvement: ", history)

        view(schur.Q, :, 1)
    elseif method==:lobpcg
        res = lobpcg(okw, okw.S, false, reshape(Pj, :, 1), 1,
                     tol=tol, C=okw.C)
        verbosity > 2 && show(io, res)

        view(res.X, :, 1)
    else
        throw(ArgumentError("Unknown diagonalization method $(method)"))
    end

    copyto!(Pj, ϕ)
end
