using QuadraticHamiltonians
using Documenter

DocMeta.setdocmeta!(QuadraticHamiltonians, :DocTestSetup, :(using QuadraticHamiltonians); recursive=true)

makedocs(;
    modules=[QuadraticHamiltonians],
    authors="cometscome <cometscome@gmail.com> and contributors",
    sitename="QuadraticHamiltonians.jl",
    format=Documenter.HTML(;
        canonical="https://cometscome.github.io/QuadraticHamiltonians.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cometscome/QuadraticHamiltonians.jl",
    devbranch="main",
)
