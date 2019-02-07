using Documenter
using SCF

makedocs(
    sitename = "SCF",
    modules = [SCF],
    pages = [
        "Home" => "index.md",
        "Quantum Systems" => "quantum_systems.md",
        "Fock Operators" => "fock_operators.md",
        "Self-Consistent Iteration" => "self_consistent_iteration.md"
    ],
    assets = ["assets/latex.js"],
)

deploydocs(repo = "github.com/JuliaAtoms/SCF.jl.git")
