using Roots
using PrettyTables

function overlap(S, u, v, tmp)
    mul!(tmp, S, v)
    u'tmp
end

function fock_matrix_element(f, S, u, v, tmp, tmp2)
    mul!(tmp, f, v)
    overlap(S, u, tmp, tmp2)
end

krylov_wrapper(kw::KrylovWrapper) = kw
krylov_wrapper(okw::OrthogonalKrylovWrapper) = okw.A

function analyze_symmetry_orbitals(fock, P::AbstractVecOrMat{T},
                                   kws, ϵ=1e-3; verbosity=0) where T
    S = fock.S

    m = length(fock.equations)
    data = zeros(T, sum(length(s)*(length(s)+1)÷2 for s in fock.symmetries), 3)
    labels = Vector{String}()
    tmp = zeros(T, size(P,1))
    tmp2 = zeros(T, size(P,1))
    ii = 1
    should_rotate = Vector{Tuple{Int,Int}}()
    for sym ∈ fock.symmetries
        for i ∈ sym
            vPi = view(P, :, i)
            A = krylov_wrapper(kws[i])
            for j ∈ sym
                j < i && continue
                vPj = view(P, :, j)
                data[ii,1] = overlap(S, vPi, vPj, tmp)
                afb = data[ii,2] = fock_matrix_element(A, S, vPi, vPj, tmp, tmp2)
                data[ii,3] = fock_matrix_element(A, S, vPj, vPi, tmp, tmp2)
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

function rotate!(P::AbstractVecOrMat, fock::Fock, A::KrylovWrapper,
                 i::Integer, j::Integer, ϵ; verbosity=0)
    S = fock.S

    a = view(P,:,i)
    b = view(P,:,j)
    tmp = similar(a)
    tmp2 = similar(b)

    fba = fock_matrix_element(A, S, b, a, tmp, tmp2)
    faa = fock_matrix_element(A, S, a, a, tmp, tmp2)
    fbb = fock_matrix_element(A, S, b, b, tmp, tmp2)

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
