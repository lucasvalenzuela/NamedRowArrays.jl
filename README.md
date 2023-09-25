# NamedRowArrays

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://lucasvalenzuela.github.io/NamedRowArrays.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://lucasvalenzuela.github.io/NamedRowArrays.jl/dev/)
[![Build Status](https://github.com/lucasvalenzuela/NamedRowArrays.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lucasvalenzuela/NamedRowArrays.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/lucasvalenzuela/NamedRowArrays.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/lucasvalenzuela/NamedRowArrays.jl)


This package for the Julia language provides an array type (the `NamedRowArray`) that knows about its row names.
This allows for indexing by name in a straightforward manner in addition to the regular indexing possibilities.

## Installation

The latest version of the package is available for Julia 1.9 and newer versions and can be installed with:
```julia
using Pkg
Pkg.add("https://github.com/lucasvalenzuela/NamedRowArrays.jl")
```

For users of the JuliaCosmoSims repository, follow
[these instructions](https://gitlab.com/juliacosmosims/JuliaCosmoSimsDocs) to install the local registry,
after which the package can be regularly installed with:
```julia
using Pkg
Pkg.add("NamedRowArrays")
```
