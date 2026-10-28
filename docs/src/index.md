```@meta
CurrentModule = NamedRowArrays
```

# NamedRowArrays

[NamedRowArrays](https://github.com/lucasvalenzuela/NamedRowArrays.jl) provides an array type
that knows the names of its rows. Rows can be selected by name in addition to the usual indexing
by position, for example the `x`, `y` and `z` rows of a matrix of particle positions.

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

A [`NamedRowArray`](@ref) wraps a vector or matrix without copying it. The row names can be given
as a `Tuple` or a `Vector` of `Symbol`s.

```jldoctest usage
julia> using NamedRowArrays

julia> pos = NamedRowArray([0.0 1.0 2.0 3.0; 0.5 0.5 0.5 0.5; 1.0 0.0 -1.0 -2.0], (:x, :y, :z))
3×4 NamedRowMatrix(::Matrix{Float64}) with eltype Float64 with row names:
 x y z
 0.0  1.0   2.0   3.0
 0.5  0.5   0.5   0.5
 1.0  0.0  -1.0  -2.0
```

A row name as the only index selects a whole row. Together with a column index, a row name or a
vector of row names selects elements, and a `Not` of row names excludes rows:

```jldoctest usage
julia> pos[:z]
4-element Vector{Float64}:
  1.0
  0.0
 -1.0
 -2.0

julia> pos[[:x, :z], 2:3]
2×2 NamedRowMatrix(::Matrix{Float64}) with eltype Float64 with row names:
 x z
 1.0   2.0
 0.0  -1.0

julia> pos[Not(:y), 1]
2-element NamedRowVector(::Vector{Float64}) with eltype Float64 with row names:
 x z
 0.0
 1.0

julia> pos[:w]
ERROR: ArgumentError: invalid row name: w
```

Row names that the array does not have are ignored in `Not`, so that the same rows can be excluded
from arrays with different row names.

Broadcasting, `map`, `hcat` and reductions that keep every row, such as `sum(pos; dims=2)`, keep
the row names.

The rows of a matrix product are the rows of the left factor. To keep the row names of `pos` when
rotating it, give the rotation matrix the same row names:

```jldoctest usage
julia> R = NamedRowArray([0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0], NamedRowArrays.rownames(pos));

julia> R * pos
3×4 NamedRowMatrix(::Matrix{Float64}) with eltype Float64 with row names:
 x y z
 -0.5  -0.5  -0.5  -0.5
  0.0   1.0   2.0   3.0
  1.0   0.0  -1.0  -2.0
```

Without row names on the left factor, the product is a plain array.

Functions that reorder rows, such as [`reverse`](@extref `Base.reverse-Tuple{AbstractVector}`),
[`sort`](@extref Base.sort), [`sortslices`](@extref Base.sortslices),
[`circshift`](@extref Base.circshift) and [`Random.shuffle`](@extref), reorder the row names as well:

```jldoctest usage
julia> reverse(pos; dims=1)
3×4 NamedRowMatrix(::Matrix{Float64}) with eltype Float64 with row names:
 z y x
 1.0  0.0  -1.0  -2.0
 0.5  0.5   0.5   0.5
 0.0  1.0   2.0   3.0
```

Their in-place versions, such as [`sort!`](@extref Base.sort!), throw an error instead.

## Reference

```@autodocs
Modules = [NamedRowArrays]
```
