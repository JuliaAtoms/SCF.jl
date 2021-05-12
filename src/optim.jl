"""
    optimize!(fock[, ::Type{Optimizer}=BFGS; kwargs...])

Solve the Hartree–Fock problem via non-linear optimization. The
algorithm, specified by `Optimizer`, has to be a first-order algorithm
that supports optimization on a manifold, the default being
`Optim.BFGS`. The computation of the gradient is accomplished by
applying the `fock` operator onto test vectors (similar to how
[`scf!`](@ref) performs Krylov iterations).

"""
function optimize!(fun!::Function, fock::Fock, ::Type{Optimizer}=BFGS;
                   opt_iters=1000, f_tol=0.0, g_tol=1e-8,
                   scf_iters=200, scf_tol=1e-3, scf_method=:lobpcg,
                   update_mixing_coefficients=true,
                   linesearch=LineSearches.HagerZhang(),
                   verbosity=2, num_printouts=typemax(Int),
                   kwargs...) where {Optimizer<:Optim.FirstOrderOptimizer}
    P = orbitals(fock.quantum_system)
    c = coefficients(fock.quantum_system)

    nc = length(c)

    H = spzeros(eltype(c), nc, nc)
    # Kinetic energy matrix
    T = spzeros(eltype(c), nc, nc)

    kws = [KrylovWrapper(hamiltonian(eq))
           for eq in fock.equations]
    f = FockProblem(fock, P, c, H, kws)
    manif = setup_manifold(f)

    trace,tolerance,_,eng,virial,flags = setup_solver_trace(
        verbosity, opt_iters, g_tol, 0, num_printouts,
        tol_header="|g|")

    trace_callback = opt_state -> begin
        # TODO: Think about moving secular problem to optimization.
        solve_secular_problem!(H, c, fock, update_mixing_coefficients)

        fun!(P, c)

        if !isnothing(trace)
            tolerance.current = opt_state.g_norm

            Etot = opt_state.value
            energy_matrix!(T, fock, :kinetic_energy)
            ET = c'T*c
            EV = Etot-ET
            eng[1].E = real(Etot)
            eng[2].E = real(ET)
            eng[3].E = real(EV)
            virial.V = real(EV/ET)

            SolverTraces.next!(trace)
        end

        false
    end

    options = Optim.Options(iterations=opt_iters,
                            f_tol=f_tol,
                            g_tol=g_tol,
                            allow_f_increases=true,
                            callback=trace_callback)

    optimizer = Optimizer(manifold=manif, linesearch=linesearch)

    t₀ = time()

    if scf_iters > 0
        @info "Performing initial SCF iterations"
        scf!(fun!, fock;
             max_iter=scf_iters, tol=scf_tol,
             method=scf_method,
             update_mixing_coefficients=update_mixing_coefficients,
             verbosity=verbosity,
             num_printouts=typemax(Int),
             kwargs...)
    end

    isnothing(trace) || print_header(trace)
    o = @time optimize(f, (w,v) -> jac!(w,v,f), copy(P),
                       optimizer, options)
    verbosity > 1 && display(o)
    copyto!(P, o.minimizer)

    rayleigh_ritz!(P, fock, kws; verbosity=verbosity)

    elapsed = time() - t₀
    verbosity > 0 && println("Finished in $(elapsed) seconds")

    analyze_symmetry_orbitals(fock, P, kws, verbosity=verbosity)

    fock
end

optimize!(fock::Fock, args...; kwargs...) =
    optimize!((_,_)->nothing, fock, args...; kwargs...)

export optimize!
