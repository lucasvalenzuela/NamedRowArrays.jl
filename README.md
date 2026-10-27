# NamedRowArrays

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://lucasvalenzuela.github.io/NamedRowArrays.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://lucasvalenzuela.github.io/NamedRowArrays.jl/dev/)
[![Build Status](https://github.com/lucasvalenzuela/NamedRowArrays.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lucasvalenzuela/NamedRowArrays.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/lucasvalenzuela/NamedRowArrays.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/lucasvalenzuela/NamedRowArrays.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)


NamedRowArrays provides an array type that knows the names of its rows. Rows can be selected by
name in addition to the usual indexing by position, for example the `x`, `y` and `z` rows of a
matrix of particle positions.

## Installation

This package is available through the public registry
[CosmoAnalysisRegistry](https://github.com/lucasvalenzuela/CosmoAnalysisRegistry). Add the registry
once per Julia installation, together with the General registry, which provides the dependencies:
```julia
using Pkg
Pkg.Registry.add("General")
Pkg.Registry.add(url="https://github.com/lucasvalenzuela/CosmoAnalysisRegistry")
```

The package requires Julia 1.12 or newer and can then be installed with:
```julia
Pkg.add("NamedRowArrays")
```

## Usage

```julia
using NamedRowArrays

pos = NamedRowArray(rand(3, 1000), (:x, :y, :z))   # wraps the matrix without copying it
pos[:z]                     # the whole row z
pos[[:x, :z], 2:3]          # a 2×2 NamedRowMatrix with rows x and z
pos[Not(:y), 1]             # a NamedRowVector with rows x and z
pos .* 2                    # broadcasting keeps the row names

R = NamedRowArray([0 -1 0; 1 0 0; 0 0 1], NamedRowArrays.rownames(pos))  # 90° about z
R * pos                     # the product has the row names of R
```
