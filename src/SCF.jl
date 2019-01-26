module SCF

using LinearAlgebra
using ContinuumArrays
import ContinuumArrays.QuasiArrays: AbstractQuasiMatrix
using ArnoldiMethod

using SolverTraces
import SolverTraces: base_exp

using UnicodeFun
using Formatting

if VERSION < v"1.1-DEV"
    isnothing(::Nothing) = true
    isnothing(::Any) = false
end

include("krylov_wrapper.jl")
include("fock.jl")
include("solver_trace.jl")
include("utils.jl")
include("self_consistent_iteration.jl")

export Fock, scf!

end # module
