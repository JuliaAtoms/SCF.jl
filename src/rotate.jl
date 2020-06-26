using Roots
using PrettyTables

fock_matrix_element(f, S, u, v, tmp) = dot(u, S, mul!(tmp, f, v))

krylov_wrapper(kw::KrylovWrapper) = kw
krylov_wrapper(okw::OrthogonalKrylovWrapper) = okw.A

function analyze_symmetry_orbitals(fock, P::AbstractVecOrMat{T},
                                   kws, ϵ=1e-3; verbosity=0) where T
    S = fock.S

    m = length(fock.equations)
    data = zeros(T, sum(length(s)*(length(s)+1)÷2 for s in fock.symmetries), 3)
    labels = Vector{String}()
    tmp = zeros(T, size(P,1))
    ii = 1
    should_rotate = Vector{Tuple{Int,Int}}()
    for sym ∈ fock.symmetries
        for i ∈ sym
            vPi = view(P, :, i)
            A = krylov_wrapper(kws[i])
            for j ∈ sym
                j < i && continue
                vPj = view(P, :, j)
                data[ii,1] = dot(vPi, S, vPj)
                afb = data[ii,2] = fock_matrix_element(A, S, vPi, vPj, tmp)
                data[ii,3] = fock_matrix_element(A, S, vPj, vPi, tmp)
                push!(labels, "$i – $j")
                ii += 1
                i ≠ j && abs(afb) > ϵ && push!(should_rotate, (i,j))
            end
        end
    end
    if verbosity > 0
        pretty_table(hcat(labels, data), ["i – j", "⟨i|j⟩", "⟨i|𝔣|j⟩", "⟨j|𝔣|i⟩"])
    end

    should_rotate
end

function rotate!(u::AbstractVector, v::AbstractVector, η, tmp::AbstractVector)
    a = inv(√(1+η^2))
    b = η*a
    copyto!(tmp, u)
    BLAS.axpy!(b, v, tmp)
    BLAS.axpy!(-b, u, v)
    copyto!(u, tmp)
    u,v
end

function rotate!(P::AbstractVecOrMat, fock::Fock, A,
                 i::Integer, j::Integer, ϵ; verbosity=0)
    S = fock.S

    a = view(P,:,i)
    b = view(P,:,j)
    tmp = similar(a)

    fba = fock_matrix_element(A, S, b, a, tmp)
    faa = fock_matrix_element(A, S, a, a, tmp)
    fbb = fock_matrix_element(A, S, b, b, tmp)

    g = (faa-fbb)
    g′ = fba
    verbosity > 0 && @show i,j,g,g′,ϵ
    abs(g′) < ϵ && return
    η̃ = g/2g′
    # TODO: Derive this properly
    η̂ = -(η̃ + √(η̃^2 + 1))
    abs(η̂) > 2 && return
    verbosity > 0 && @show i,j,η̂

    rotate!(a, b, η̂, tmp)
end

rotate!(P::AbstractVecOrMat, fock::Fock, okws::Vector,
        i::Integer, j::Integer, ϵ;
        kwargs...) =
            rotate!(P, fock, krylov_wrapper(okws[i]),
                    i, j, ϵ; kwargs...)

function rayleigh_ritz!(P, fock, kws; verbosity=0)
    verbosity > 0 && @info "Rayleigh–Ritz rotation"
    T = eltype(P)
    S = fock.S
    tmp = zeros(T, size(P,1))
    for sym ∈ fock.symmetries
        ns = length(sym)
        if verbosity > 2
            println(repeat("-", 100))
            println("Symmetry: ", sym)
        end
        f = zeros(T, ns, ns)
        Ps = P[:,sym]
        for (ii,i) ∈ enumerate(sym)
            vPi = view(Ps, :, ii)
            A = krylov_wrapper(kws[i])
            for (jj,j) ∈ enumerate(sym)
                j < i && continue
                vPj = view(Ps, :, jj)
                f[ii,jj] = fock_matrix_element(A, S, vPi, vPj, tmp)
                f[jj,ii] = fock_matrix_element(A, S, vPj, vPi, tmp)
            end
        end
        ee = eigen(f)
        mul!(view(P,:,sym)', ee.vectors', Ps')
        if verbosity > 2
            println("Fock submatrix:")
            display(f)
            println("Energies:")
            display(ee.values')
            display(ee.vectors)
            println()
        end
    end
    update!(fock)
    true
end

function rotate_orbitals!(P, fock, kws;
                          verbosity=0,
                          rotate_orbitals=true, rotϵ=1e-3,
                          rotation_method=:pairwise,
                          kwargs...)
    did_rotate = if rotate_orbitals
        if rotation_method == :pairwise
            should_rotate=analyze_symmetry_orbitals(fock, P, kws, rotϵ,
                                                    verbosity=verbosity)
            for (i,j) in should_rotate
                rotate!(P, fock, kws, i, j, rotϵ, verbosity=verbosity)
            end
            !isempty(should_rotate)
        elseif rotation_method == :eigen
            rayleigh_ritz!(P, fock, kws; verbosity=verbosity)
        else
            throw(ArgumentError("Unknown orbital rotation method $(rotation_method)"))
        end
    else
        false
    end
    did_rotate && verbosity > 2 &&
        analyze_symmetry_orbitals(fock, P, kws, verbosity=verbosity)
    did_rotate
end
