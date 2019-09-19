mutable struct FockProblem{F,O,C,M,KWS}
    fock::F
    P::O
    c::C
    H::M
    kws::KWS
end

function set!(f::FockProblem, v)
    copyto!(f.P, v)
    update!(f.fock)
end

function (f::FockProblem)()
    energy_matrix!(f.H, f.fock)
    f.c'f.H*f.c
end

function (f::FockProblem)(v)
    set!(f, v)
    f()
end

function jac!(w, v, f::FockProblem)
    set!(f, v)
    # Threads.@threads
    for j = 1:length(f.fock.equations)
        mul!(view(w,:,j), f.kws[j], view(v,:,j))
    end
end

"""
    optimize!(fock[, ::Type{Optimizer}=BFGS; kwargs...])

Solve the Hartree–Fock problem via non-linear optimization. The
algorithm, specified by `Optimizer`, has to be a first-order algorithm
that supports optimization on a manifold, the default being
`Optim.BFGS`. The computation of the gradient is accomplished by
applying the `fock` operator onto test vectors (similar to how
[`scf!`](@ref) performs Krylov iterations).

"""
function optimize!(fock::Fock, ::Type{Optimizer}=BFGS;
                   opt_iters=1000, g_tol=1e-8,
                   scf_iters=0,
                   verbosity=2, num_printouts=typemax(Int),
                   kwargs...) where {Optimizer<:Optim.FirstOrderOptimizer}
    P = orbitals(fock.quantum_system)
    c = coefficients(fock.quantum_system)

    m,n = size(P)
    nc = length(c)

    H = spzeros(nc, nc)
    # Kinetic energy matrix
    T = spzeros(nc, nc)

    if scf_iters > 0
        @info "Performing initial SCF iterations"
        scf!(fock;
             max_iter=scf_iters, verbosity=2,
             num_printouts=typemax(Int),
             kwargs...)
    end

    trace,tolerance,_,eng,virial = setup_solver_trace(
        verbosity, opt_iters, g_tol, 0, num_printouts)

    trace_callback = opt_state -> begin
        # TODO: Think about moving secular problem to optimization.
        solve_secular_problem!(H, c, fock)

        if !isnothing(trace)
            tolerance.current = opt_state.g_norm

            Etot = opt_state.value
            energy_matrix!(T, fock, :kinetic_energy)
            ET = c'T*c
            EV = Etot-ET
            eng[1].E = Etot
            eng[2].E = ET
            eng[3].E = EV
            virial.V = EV/ET

            SolverTraces.next!(trace)
        end

        false
    end

    options = Optim.Options(iterations=opt_iters, g_tol=g_tol,
                            allow_f_increases=true,
                            callback=trace_callback)

    manif = setup_manifold(fock, m, n)

    kws = [KrylovWrapper(hamiltonian(eq))
           for eq in fock.equations]
    f = FockProblem(fock, P, c, H, kws)

    optimizer = Optimizer(manifold=manif)

    print_header(trace)
    o = @time optimize(f, (w,v) -> jac!(w,v,f), copy(P),
                       optimizer, options)
    verbosity > 1 && display(o)
    copyto!(P, o.minimizer)

    fock
end

export optimize!
