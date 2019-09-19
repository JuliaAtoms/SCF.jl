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

function analyze_symmetry_orbitals(fock, P::AbstractVecOrMat{T},
                                   okws; verbosity=0) where T
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
            A = okws[i].A
            for j ∈ sym
                j < i && continue
                vPj = view(P, :, j)
                data[ii,1] = overlap(S, vPi, vPj, tmp)
                afb = data[ii,2] = fock_matrix_element(A, S, vPi, vPj, tmp, tmp2)
                data[ii,3] = fock_matrix_element(A, S, vPj, vPi, tmp, tmp2)
                push!(labels, "$i – $j")
                ii += 1
                i ≠ j && abs(afb) > 1e-3 && push!(should_rotate, (i,j))
            end
        end
    end
    if verbosity > 3
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

function rotate!(P::AbstractVecOrMat, fock::Fock, okws, i::Integer, j::Integer)
    S = fock.S
    A = okws[i].A

    u = view(P,:,i)
    v = view(P,:,j)
    a = similar(u)
    b = similar(u)
    tmp = similar(u)
    tmp2 = similar(u)

    f = η -> begin
        copyto!(a, u)
        copyto!(b, v)
        rotate!(a, b, η, tmp)
        fock_matrix_element(A, S, a, b, tmp, tmp2)
    end
    ηlo,ηhi = (-1,1)
    flo,fhi = (f(ηlo),f(ηhi))
    η = if sign(flo) ≠ sign(fhi)
        find_zero(f, (ηlo,ηhi), Bisection())
    else
        f₀ = f(0)
        if abs(flo) < f₀
            ηlo
        elseif abs(fhi) < f₀
            ηhi
        else
            0
        end
    end
    rotate!(u, v, η, tmp)
end
