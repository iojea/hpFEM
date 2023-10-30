using FEMhp
using Documenter

DocMeta.setdocmeta!(FEMhp, :DocTestSetup, :(using FEMhp); recursive=true)

makedocs(;
    modules=[FEMhp],
    authors="Ignacio Ojea",
    repo="https://github.com/iojea/FEMhp.jl/blob/{commit}{path}#{line}",
    sitename="FEMhp.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://iojea.github.io/FEMhp.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/iojea/FEMhp.jl",
    devbranch="main",
)
