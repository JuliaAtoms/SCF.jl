mutable struct KrylovWrapper{T,B<:AbstractQuasiMatrix,Hamiltonian}
    R::B
    hamiltonian::Hamiltonian
end

Base.eltype(A::KrylovWrapper{T}) where T = T
Base.size(A::KrylovWrapper) = (size(A.R,2),size(A.R,2))
Base.size(A::KrylovWrapper, i) = size(A)[i]

Base.show(io::IO, kw::KrylovWrapper{T,B,Hamiltonian}) where {T,B,Hamiltonian} =
    write(io, "KrylovWrapper{$T} of size $(size(kw))")
