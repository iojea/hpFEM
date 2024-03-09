using hpFEM
using Documenter

DocMeta.setdocmeta!(hpFEM, :DocTestSetup, :(using hpFEM); recursive=true)

makedocs(;
    modules=[hpFEM],
    authors="Ignacio Ojea",
    repo="https://github.com/iojea/hpFEM.jl/blob/{commit}{path}#{line}",
    sitename="FEMhp.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://iojea.github.io/hpFEM.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/iojea/hpFEM.jl",
    devbranch="main",
)
