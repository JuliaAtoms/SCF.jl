using Documenter
using SCF

makedocs(
    sitename = "SCF",
    modules = [SCF],
    assets = ["assets/latex.js"],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
