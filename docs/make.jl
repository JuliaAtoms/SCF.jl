using SCF
using Documenter

makedocs(;
    modules=[SCF],
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com>",
    repo="https://github.com/JuliaAtoms/SCF.jl/blob/{commit}{path}#{line}",
    sitename="SCF.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="juliaatoms.org",
        assets=String["assets/latex.js"],
        mathengine = Documenter.MathJax()
    ),
    pages=[
        "Home" => "index.md",
        "Quantum Systems" => "quantum_systems.md",
        "Fock Operators" => "fock_operators.md",
        "Self-Consistent Iteration" => "self_consistent_iteration.md"
    ],
    doctest=false,
)

deploydocs(;
    repo="github.com/JuliaAtoms/SCF.jl",
)
