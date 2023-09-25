module NamedRowArrays

using Base: @propagate_inbounds

export NamedRowArray, NamedRowVector, NamedRowMatrix

"""
A NamedRowArray is an AbstractArray that wraps another AbstractArray (only
AbstractVector and AbstractMatrix, however) and adds row names to the first array
dimension. AxisArrays can be indexed by using the named rows as an alternative to
positional indexing by dimension.

### Type parameters

The NamedRowArray contains several type parameters:

```julia
struct AxisArray{T,N,D,R} <: AbstractArray{T,N}
```
- `T` : the elemental type of the AbstractArray
- `N` : the number of dimensions
- `D` : the type of the wrapped AbstractArray
- `R` : the names of the rows, as a NTuple{N, Symbol}

### Constructors

When constructing a NamedRowArray, note that no copy of the array is made.

```julia
NamedRowArray(A::AbstractArray, rownames::Tuple)
NamedRowArray(A::AbstractArray, rownames::AbstractVector)
```

### Indexing

NamedRowArrays can be indexed in the same way as normal AbstractArrays.
In addition, one can also use the rownames instead of the numerical indices.
Since only the rows can have names, it is also possible to obtain one or more
rows of a Matrix by only passing one or more row names as the index.

### Examples

For a Vector:
```julia
v = NamedRowArray(rand(3), [:a, :b, :c])
v[1] # scalar value
v[:a] # scalar value
v[2:3] # 2-element NamedRowVector
v[[:a, :b]] # 2-element NamedRowVector
```

For a Matrix:
```julia
a = NamedRowArray(rand(3, 10), [:a, :b, :c])
a[1, 3] # scalar value
a[:a, 3] # scalar value
a[:a, 2:3] # 2-element Vector
a[:a] # whole row: 10-element Vector
a[[:a, :b]] # two rows: 2×10 NamedRowMatrix
a[:, 3] # one column: 3-element NamedRowVector
a[:, 2:3] # two columns: 3×2 NamedRowMatrix
```
"""
struct NamedRowArray{T,N,D,R} <: AbstractArray{T,N}
    data::D # D <:AbstractArray, enforced in constructor to avoid dispatch bugs
    rownames::R # R<:NTuple{N, Symbol}
    NamedRowArray{T,N,D,R}(data::AbstractArray{T,N}, rownames::NTuple) where {T,N,D,R} = new{T,N,D,R}(data, rownames)
end

"""
    NamedRowMatrix{T,D,R}

Alias for [`NamedRowArray{T,2,D,R}`](@ref NamedRowArray).
"""
const NamedRowMatrix{T,D,R} = NamedRowArray{T,2,D,R}

"""
    NamedRowVector{T,D,R}

Alias for [`NamedRowArray{T,1,D,R}`](@ref NamedRowArray).
"""
const NamedRowVector{T,D,R} = NamedRowArray{T,1,D,R}

NamedRowArray(A::AbstractArray, rownames::AbstractArray) = NamedRowArray(A, Tuple(rownames))
NamedRowVector(A::AbstractVector, rownames) = NamedRowArray(A, Tuple(rownames))
NamedRowMatrix(A::AbstractMatrix, rownames) = NamedRowArray(A, Tuple(rownames))


function NamedRowArray(A::D, rownames::R) where {T,N,D<:AbstractArray{T,N},R<:NTuple}
    A isa AbstractVecOrMat || throw(ArgumentError("the array has to be an AbstractVecOrMat"))
    eltype(rownames) === Symbol || throw(ArgumentError("the rownames have to be of type Symbol"))
    length(rownames) == _size(A) || throw(ArgumentError("the length of the rownames must match the corresponding row length of the data"))
    return NamedRowArray{T,N,D,R}(A, rownames)
end

_size(A::AbstractVector) = length(A)
_size(A::AbstractMatrix) = size(A, 1)

# Base definitions that aren't provided by AbstractArray
@inline Base.size(A::NamedRowArray) = size(A.data)
@inline Base.axes(A::NamedRowArray) = axes(A.data)
Base.convert(::Type{Array{T,N}}, A::NamedRowArray{T,N}) where {T,N} = convert(Array{T,N}, A.data)
Base.parent(A::NamedRowArray) = A.data
Base.similar(A::NamedRowArray, ::Type{S}) where {S} = (d = similar(A.data, S); NamedRowArray(d, A.rownames))
Base.similar(A::NamedRowArray, ::Type{S}, dims::Dims{N}) where {S,N} = similar(A.data, S, dims)

# A simple display method to include axis information. It might be nice to
# eventually display the axis labels alongside the data array, but that is
# much more difficult.
function Base.summary(io::IO, A::NamedRowArray)
    _summary(io, A)
    print(io, join(':' .* String.(A.rownames), ", "))
    println(io)
    print(io, "for the ", summary(A.data))
end
_summary(io, _::NamedRowVector{T,N}) where {T,N} = println(io, "NamedRowVector{$T,...} with row names:")
_summary(io, _::NamedRowMatrix{T,N}) where {T,N} = println(io, "NamedRowMatrix{$T,...} with row names:")



# Simple scalar indexing where we just set or return scalars
@propagate_inbounds Base.getindex(A::NamedRowArray, idxs::Int...) = A.data[idxs...]
@propagate_inbounds Base.setindex!(A::NamedRowArray, v, idxs::Int...) = (A.data[idxs...] = v)

const Idx = Union{Colon,AbstractArray{Int},BitVector,AbstractArray{Bool}}


reaxis(A::NamedRowArray, I::Idx) = A.rownames[I]
reaxis(A::NamedRowMatrix, I::Idx, _) = A.rownames[I]

@propagate_inbounds function Base.getindex(A::NamedRowArray, idxs::Idx...)
    NamedRowArray(A.data[idxs...], reaxis(A, idxs...))
end

@propagate_inbounds Base.getindex(A::NamedRowMatrix, I::Idx) = A.data[I]
@propagate_inbounds Base.getindex(A::NamedRowMatrix, idx::Int, I::Idx) = A.data[idx, I]
@propagate_inbounds Base.getindex(A::NamedRowMatrix, idx::Int, I::Int) = A.data[idx, I]
@propagate_inbounds Base.getindex(A::NamedRowMatrix, idx::Idx, I::Int) = NamedRowVector(A.data[idx, I], reaxis(A, idx))

# To resolve ambiguities, we need several definitions
@propagate_inbounds Base.view(A::NamedRowArray, idxs::Idx...) = NamedRowArray(view(A.data, idxs...), reaxis(A, idxs...))
@propagate_inbounds Base.view(A::NamedRowMatrix, idx::Int, I::Idx) = view(A.data, idx, I)

# Setindex is so much simpler. Just assign it to the data:
@propagate_inbounds Base.setindex!(A::NamedRowArray, v, idxs::Union{Idx,Int}...) = (A.data[idxs...] = v)


# Indexing by column names
@propagate_inbounds Base.getindex(A::NamedRowVector, idx::Symbol) = getindex(A, to_index(A, idx))
@propagate_inbounds Base.getindex(A::NamedRowMatrix, idx::Symbol, I=(:)) = getindex(A, to_index(A, idx), I)
@propagate_inbounds Base.getindex(A::NamedRowVector, idxs::AbstractVector{Symbol}) = getindex(A, to_index.((A,), idxs))
@propagate_inbounds Base.getindex(A::NamedRowMatrix, idxs::AbstractVector{Symbol}, I=(:)) = getindex(A, to_index.((A,), idxs), I)

@propagate_inbounds Base.setindex!(A::NamedRowVector, v, idx::Symbol) = setindex!(A, v, to_index(A, idx))
@propagate_inbounds Base.setindex!(A::NamedRowMatrix, v, idx::Symbol, I=(:)) = setindex!(A, v, to_index(A, idx), I)
@propagate_inbounds Base.setindex!(A::NamedRowVector, v, idxs::AbstractVector{Symbol}) = setindex!(A, v, to_index.((A,), idxs))
@propagate_inbounds Base.setindex!(A::NamedRowMatrix, v, idxs::AbstractVector{Symbol}, I=(:)) = setindex!(A, v, to_index.((A,), idxs), I)


function to_index(A::NamedRowArray, idx::Symbol)
    ind = findfirst(==(idx), A.rownames)
    isnothing(ind) && throw(ArgumentError("invalid index: $idx"))
    return ind
end

# Cartesian iteration
Base.eachindex(A::NamedRowArray) = eachindex(A.data)

end # module
