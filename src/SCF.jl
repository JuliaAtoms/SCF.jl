module SCF

using LinearAlgebra
using ArnoldiMethod

using SolverTraces
import SolverTraces: base_exp

using UnicodeFun
using Formatting

if VERSION < v"1.1-DEV"
    isnothing(::Nothing) = true
    isnothing(::Any) = false
end

include("quantum_systems.jl")
include("krylov_wrapper.jl")
include("fock.jl")
include("solver_trace.jl")
include("utils.jl")
include("self_consistent_iteration.jl")

export Fock, scf!

end # module
