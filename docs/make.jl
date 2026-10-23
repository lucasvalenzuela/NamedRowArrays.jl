using NamedRowArrays
using Documenter
using DocumenterInterLinks

DocMeta.setdocmeta!(NamedRowArrays, :DocTestSetup, :(using NamedRowArrays); recursive=true)

links = InterLinks("Julia" => "https://docs.julialang.org/en/v1/")

makedocs(;
    plugins=[links],
    modules=[NamedRowArrays],
    authors="Lucas Valenzuela",
    sitename="NamedRowArrays.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://lucasvalenzuela.github.io/NamedRowArrays.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/lucasvalenzuela/NamedRowArrays.jl",
    devbranch="main",
)
