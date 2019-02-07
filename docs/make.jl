using Documenter
using SCF

makedocs(
    sitename = "SCF",
    modules = [SCF],
    assets = ["assets/latex.js"],
)

deploydocs(repo = "github.com/JuliaAtoms/SCF.jl.git")
