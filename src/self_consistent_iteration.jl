"""
    circshiftpush!(v, e)

Shifts all elements of the vector `v` one step to the left and sets
`e` as the last element.
"""
function circshiftpush!(v::V, e::T) where {T,V<:AbstractVector{T}}
    for i = 1:length(v)-1
        v[i] = v[i+1]
    end
    v[end] = e
    v
end

"""
    ismonotonous(v[, comp=isless])

Returns true if the magnitude of every value fulfils `comp(cur,prev)`.
"""
function ismonotonous(v::V, comp=isless) where {T,V<:AbstractVector{T}}
    for i = 2:length(v)
        comp(abs(v[i]), abs(v[i-1])) || return false
    end
    true
end


"""
    scf!([fun!, ]fock[; ω=0, max_iter=200, tol=1e-8, verbosity=1])

Perform the SCF iterations for the `fock` operator. One iteration is
performed by [`scf_iteration!`](@ref). The optional `fun!(P̃,c̃)`
argument allows for extra operations to be performed every SCF cycle
(such as plotting, etc).

`ω` is a relaxation parameter; the orbitals and mixing coefficients
are updated as `wᵢ₊₁ = (1-ω)w̃ + ωwᵢ` where `w̃` is the solution in the
current iteration and `wᵢ` the previous solution. The default (`ω=0`)
is to only use this solution. A larger value of `ω` makes it easier to
achieve convergence in oscillatory problems, but a too large value
decreases the convergence rate. It is therefore desireable to reduce
the value of `ω` if the convergence is deemed to be monotonous, the
criterion for which is whether the change over the `monotonous_window`
last iterations is steadily decreasing. If this is the case, `ω` will
be reduced by multiplying by `ωfactor`. Conversely, if the change is
non-monotonous, `ω` will be increased, but not beyond `ωmax`.

The SCF procedure continues until either the amount of iterations
equals `max_iter` or the change in the coefficients is below `tol`.
"""
function scf!(fun!::Function, fock::Fock{Q};
              max_iter=200, tol=1e-8,
              ω=0, monotonous_window=10, ωfactor=0.9, ωmax=0.999,
              double_counted_total_energy=false,
              verbosity=1, num_printouts=min(max_iter,10),
              kwargs...) where Q
    trace,tolerance,relaxation,eng,virial =
        setup_solver_trace(verbosity, max_iter, tol, ω, num_printouts)

    # Views of the orbitals and mixing coefficients in the fock object.
    P = orbitals(fock.quantum_system)
    c = coefficients(fock.quantum_system)

    # Current estimates
    P̃ = copy(P)
    c̃ = copy(c)

    norb = size(P,2)
    nc = length(c)

    ΔP = [Inf for j in 1:norb]
    Δc = [Inf for j in 1:nc]

    Δhistory = zeros(monotonous_window)

    if verbosity ≥ 2
        nc > 1 && print("Multi-Configurational ")
        println("Self-Consistent-Field calculation of")
        print("- ")
        display(fock.quantum_system)
        println("- Maximum amount of iterations: $(max_iter)")
        tb,te = base_exp(tol)
        println("- Stopping tolerance: ", format(tolerance.tol_fmt, tb, to_superscript(te)))
        ω != 0 && println("- Successive relaxation: ω = $(ω)")
        println()
    end

    !isnothing(trace) && print_header(trace)
    t₀ = time()

    # Energy matrix
    H = spzeros(nc, nc)
    solve_secular_problem!(H, c, fock)
    # Energy matrix for double-counted terms
    Hd = spzeros(nc, nc)
    # Kinetic energy matrix
    T = spzeros(nc, nc)

    if verbosity > 2
        println("Initial energy matrix")
        display(Matrix(H))
        println()
        println("Initial mixing coefficients:")
        display(c)
        println()
    end

    # Prepare all orbital integrals before the first iteration.
    update!(fock.equations; verbosity=max(0,verbosity-4))

    for i = 1:max_iter
        scf_iteration!(fock, P̃, H, c̃;
                       update_mixing_coefficients=i>1,
                       verbosity=verbosity-2, kwargs...)
        fun!(P̃, c̃)

        for j = 1:norb
            ΔP[j] = 1.0 - dot(view(P̃,:,j),P[:,j])/norm(P[:,j])^2
        end
        Δc .= c̃ .- c

        aΔ = norm(ΔP) + norm(Δc)
        isnothing(tolerance) || (tolerance.current = aΔ)

        if ω == 0
            copyto!(P, P̃)
            copyto!(c, c̃)
        else
            # Relaxation; see e.g. p. 490 of
            #
            # - Fischer, C. F., & Guo, W. (1990). Spline Algorithms for the
            #   Hartree-Fock Equation for the Helium Ground State. Journal of
            #   Computational Physics, 90(2),
            #   486–496. http://dx.doi.org/10.1016/0021-9991(90)90176-2
            lmul!(ω, P)
            lmul!(1-ω, P̃)
            P[:] += P̃[:] # Is this efficient?
            lmul!(ω, c)
            lmul!(1-ω, c̃)
            c[:] += c̃[:] # Is this efficient?
            normalize!(c)
        end

        circshiftpush!(Δhistory, aΔ)
        if ismonotonous(Δhistory)
            ω *= ωfactor
        else
            ω /= ωfactor
        end
        ω = clamp(ω, 0.0, ωmax)

        isnothing(relaxation) || (relaxation.ω = ω)
        if !isnothing(eng)
            Etot = if double_counted_total_energy
                energy_matrix!(Hd, fock, :double_counted_energy)
                sum(energy.(fock.equations)) - c'Hd*c
            else
                c'H*c
            end
            energy_matrix!(T, fock, :kinetic_energy)
            ET = c'T*c
            EV = Etot-ET
            eng[1].E = Etot
            eng[2].E = ET
            eng[3].E = EV
            virial.V = EV/ET
        end
        SolverTraces.next!(trace)

        if aΔ < tol
            println()
            break
        end
    end
    elapsed = time() - t₀
    verbosity > 0 && println("Finished in $(elapsed) seconds")

    norm(ΔP) + norm(Δc) > tol && @warn "Desired tolerance $(tol) not reached in $(max_iter) iterations"

    fock
end

scf!(fock::Fock{Q}; kwargs...) where Q =
    scf!((_,_)->nothing, fock; kwargs...)


"""
    scf_iteration!(fock, P, c[; verbosity=0, method=:arnoldi, σ=nothing, tol=1e-10,
                   update_mixing_coefficients=true])

Perform one step of the SCF iteration. This will

1) Improve each of the orbitals in turn, using the values of the
   orbitals `P` and mixing coefficients `c` from the previous step as
   input.

2) Solve the secular problem to find new values for the mixing
   coefficients, `c`, unless `update_mixing_coefficients==false`.

"""
function scf_iteration!(fock::Fock{Q,E}, P::M, H::HM, c::V;
                        verbosity=0, method=:arnoldi,
                        tol=1e-10,
                        update_mixing_coefficients=true,
                        kwargs...) where {Q,E,M<:AbstractVecOrMat,
                                          HM<:AbstractMatrix,V<:AbstractVector}
    verbosity > 0 && println("Improving orbitals using $(method)")

    for (j,eq) in enumerate(fock.equations)
        print_block() do io
            verbosity > 1 && println(io, eq)

            vPj = view(P,:,j)

            solve_orbital_equation!(vPj, eq, method, tol;
                                    io=io, verbosity=verbosity,
                                    kwargs...)

            # Normalize eigenvector and rotate it such that the
            # largest lobe is positive.
            norm_rot!(fock, vPj)

            verbosity > 2 &&
                println(io, "Change in equation $j: ",
                        norm(vPj - eq.ϕ.args[2]))
        end
    end

    # energy_matrix! requires fresh integrals to be calculated
    # correctly, hence we update the integrals after improving the
    # orbitals.
    update!(fock.equations; verbosity=max(0,verbosity-2))

    update_mixing_coefficients &&
        solve_secular_problem!(H, c, fock, tol=tol)

    fock
end
