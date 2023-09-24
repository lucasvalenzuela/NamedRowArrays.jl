using NamedRowArrays
using Documenter

DocMeta.setdocmeta!(NamedRowArrays, :DocTestSetup, :(using NamedRowArrays); recursive=true)

makedocs(;
    modules=[NamedRowArrays],
    authors="Lucas Valenzuela",
    repo="https://github.com/lucasvalenzuela/NamedRowArrays.jl/blob/{commit}{path}#{line}",
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
