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
