"""
    gershgorin_disc(i, row)

Return the [Gershgorin
disc](https://en.wikipedia.org/wiki/Gershgorin_circle_theorem) for the
`i`th `row` of a matrix.
"""
function gershgorin_disc(i, row::AbstractVector)
    aᵢᵢ = row[i]
    aᵢᵢ, (sum(abs, row) - abs(aᵢᵢ))
end

"""
    gershgorin_discs(A)

Compute all [Gershgorin
discs](https://en.wikipedia.org/wiki/Gershgorin_circle_theorem) of a
matrix `A`.
"""
gershgorin_discs(A::AbstractMatrix) =
    map(((i,row),) -> gershgorin_disc(i,row), enumerate(eachrow(A)))

"""
    gershgorin_bounds(A)

For a Hermitian matrix `A`, compute the bounds on the real axis of the
eigenspectrum. NB that this bound will only be useful for diagonally
dominant `A`s.
"""
function gershgorin_bounds(A::AbstractMatrix)
    ishermitian(A) ||
        throw(ArgumentError("Not implemented for non-Hermitian matrices"))
    bounds = Inf..(-Inf)
    for (i,row) in enumerate(eachrow(A))
        a,r = gershgorin_disc(i,row)
        bounds = bounds ∪ ((a-r)..(a+r))
    end
    bounds
end
