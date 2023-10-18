module NamedRowArrays

using Base: @propagate_inbounds, to_index
using InvertedIndices

export NamedRowArray, NamedRowVector, NamedRowMatrix, Not

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

### Inverted Indexing

Finally, it is also possible to use inverted indexes with NamedRowArrays by using the exported
`Not` struct from InvertedIndices.jl. In contrast to the usual behavior of `Not`, using row
names that do not exist will not throw an error, but will silently ignore them.

```julia
v[Not(:a)] # 3-element NamedRowVector
v[Not(:a, :c)] # 2-element NamedRowVector

a[Not(:a)] # three rows: 3×10 NamedRowMatrix
a[Not(:a), 2:3] # three rows, two columns: 3×2 NamedRowMatrix
a[Not(:a), 2] # 3-element NamedRowVector
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
Base.reinterpret(::Type{T}, A::NamedRowArray) where {T} = NamedRowArray(reinterpret(T, A.data), A.rownames)


# A simple display method to include axis information. It might be nice to
# eventually display the axis labels alongside the data array, but that is
# much more difficult.
function Base.summary(io::IO, A::NamedRowArray)
    _summary(io, A)
    # print(io, join(':' .* String.(A.rownames), ", "))
    _summary(io, A.rownames)
    println(io)
    dstmp = displaysize(io)
    ds = (dstmp[1] - 5, dstmp[2] - 5)
    # print(io, "for the ", summary(IOContext(io, :displaysize => ds), A.data))
    print(IOContext(io, :displaysize => ds), "for the ", summary(A.data))
end
_summary(io, _::NamedRowVector{T,N}) where {T,N} = println(io, "NamedRowVector{$T,...} with row names:")
_summary(io, _::NamedRowMatrix{T,N}) where {T,N} = println(io, "NamedRowMatrix{$T,...} with row names:")
function _summary(io, vals::NTuple)
    width = displaysize(io)[2]

    str = join(String.(vals), " ")
    if width - 1 < length(str)
        ind = findprev(' ', str, width - 2)
        print(io, " ", @view(str[1:min(end,ind)]), "…")
    else
        print(io, " ", str)
    end
end

# adapted from arrayshow.jl in Base Julia
function Base.show(io::IO, mime::MIME"text/plain", A::NamedRowArray)
    if isempty(A) && get(io, :compact, false)::Bool
        return show(io, A)
    end
    # 0) show summary before setting :compact
    summary(io, A)
    isempty(A) && return
    print(io, ":")
    Base.show_circular(io, A) && return

    # 1) compute new IOContext
    if !haskey(io, :compact) && length(axes(A, 2)) > 1
        io = IOContext(io, :compact => true)
    end
    if get(io, :limit, false)::Bool && eltype(A) === Method
        # override usual show method for Vector{Method}: don't abbreviate long lists
        io = IOContext(io, :limit => false)
    end

    dstmp = displaysize(io)
    if get(io, :limit, false)::Bool && dstmp[1]-6 <= 0
        return print(io, " …")
    else
        println(io)
    end

    # 2) update typeinfo
    #
    # it must come after printing the summary, which can exploit :typeinfo itself
    # (e.g. views)
    # we assume this function is always called from top-level, i.e. that it's not nested
    # within another "show" method; hence we always print the summary, without
    # checking for current :typeinfo (this could be changed in the future)
    ds = (dstmp[1]-2, dstmp[2])
    io = IOContext(io, :typeinfo => eltype(A), :displaysize => ds)

    # 2) show actual content
    recur_io = IOContext(io, :SHOWN_SET => A)
    Base.print_array(recur_io, A)
end



# Simple scalar indexing where we just set or return scalars
@propagate_inbounds Base.getindex(A::NamedRowArray, idxs::Int...) = A.data[idxs...]
@propagate_inbounds Base.setindex!(A::NamedRowArray, v, idxs::Int...) = (A.data[idxs...] = v)

const Idx = Union{Colon,AbstractArray{Int},BitVector,AbstractArray{Bool}}


reaxis(A::NamedRowArray, I::Idx) = A.rownames[I]
reaxis(A::NamedRowMatrix, I::Idx, _) = A.rownames[I]
reaxis(A::NamedRowArray, I::InvertedIndex) = collect(A.rownames)[I]
reaxis(A::NamedRowMatrix, I::InvertedIndex, _) = collect(A.rownames)[I]

@propagate_inbounds function Base.getindex(A::NamedRowArray, idxs::Idx...)
    NamedRowArray(A.data[idxs...], reaxis(A, idxs...))
end

@propagate_inbounds Base.getindex(A::NamedRowMatrix, I::Idx) = A.data[I]
@propagate_inbounds Base.getindex(A::NamedRowMatrix, idx::Int, I::Idx) = A.data[idx, I]
@propagate_inbounds Base.getindex(A::NamedRowMatrix, idx::Int, I::Int) = A.data[idx, I]
@propagate_inbounds Base.getindex(A::NamedRowMatrix, idx::Idx, I::Int) = NamedRowVector(A.data[idx, I], reaxis(A, idx))

@propagate_inbounds Base.view(A::NamedRowVector, idx::Symbol) = view(A.data, to_index(A, idx))
@propagate_inbounds Base.view(A::NamedRowMatrix, idx::Symbol, I=(:)) = view(A.data, to_index(A, idx), I)
@propagate_inbounds Base.view(A::NamedRowVector, idxs::AbstractVector{Symbol}) = view(A, to_index.((A,), idxs))
@propagate_inbounds Base.view(A::NamedRowMatrix, idxs::AbstractVector{Symbol}, I=(:)) = view(A, to_index.((A,), idxs), I)

# To resolve ambiguities, we need several definitions
@propagate_inbounds Base.view(A::NamedRowArray, idxs::Idx...) = NamedRowArray(view(A.data, idxs...), reaxis(A, idxs...))
@propagate_inbounds Base.view(A::NamedRowVector, idx::Int) = view(A.data, idx)
@propagate_inbounds Base.view(A::NamedRowMatrix, idx::Int, I) = view(A.data, idx, I)
@propagate_inbounds Base.view(A::NamedRowMatrix, idx::Idx, I::Int) = NamedRowArray(view(A.data, idx, I), reaxis(A, idx))

# Setindex is so much simpler. Just assign it to the data:
@propagate_inbounds Base.setindex!(A::NamedRowArray, v, idxs::Union{Idx,Int}...) = (A.data[idxs...] = v)


# Indexing by column names
@propagate_inbounds Base.getindex(A::NamedRowVector, idx::Symbol) = getindex(A, to_index(A, idx))
@propagate_inbounds Base.getindex(A::NamedRowMatrix, idx::Symbol, I=(:)) = getindex(A, to_index(A, idx), I)
@propagate_inbounds Base.getindex(A::NamedRowVector, idxs::AbstractVector{Symbol}) = getindex(A, to_index.((A,), idxs))
@propagate_inbounds Base.getindex(A::NamedRowMatrix, idxs::AbstractVector{Symbol}, I=(:)) = getindex(A, to_index.((A,), idxs), I)
@inline Base.getindex(A::NamedRowMatrix, idx, idx2::Symbol) = throw(ArgumentError("only named rows are allowed, not columns"))
@inline Base.getindex(A::NamedRowMatrix, idx::Symbol, idx2::Symbol) = throw(ArgumentError("only named rows are allowed, not columns"))

@propagate_inbounds Base.setindex!(A::NamedRowVector, v, idx::Symbol) = setindex!(A, v, to_index(A, idx))
@propagate_inbounds Base.setindex!(A::NamedRowMatrix, v, idx::Symbol, I=(:)) = setindex!(A, v, to_index(A, idx), I)
@propagate_inbounds Base.setindex!(A::NamedRowVector, v, idxs::AbstractVector{Symbol}) = setindex!(A, v, to_index.((A,), idxs))
@propagate_inbounds Base.setindex!(A::NamedRowMatrix, v, idxs::AbstractVector{Symbol}, I=(:)) = setindex!(A, v, to_index.((A,), idxs), I)
@propagate_inbounds @inline Base.setindex!(A::NamedRowMatrix, v, idx, idx2::Symbol) = throw(ArgumentError("only named rows are allowed, not columns"))
@propagate_inbounds @inline Base.setindex!(A::NamedRowMatrix, v, idx::Symbol, idx2::Symbol) = throw(ArgumentError("only named rows are allowed, not columns"))


# inverted indices
@inline Base.getindex(A::NamedRowVector, I::InvertedIndex{<:AbstractVector{Symbol}}) = getindex(A, Not(I.skip...))
@inline Base.getindex(A::NamedRowMatrix, I::InvertedIndex{<:AbstractVector{Symbol}}) = getindex(A, Not(I.skip...))
@inline Base.getindex(A::NamedRowMatrix, I::InvertedIndex{<:AbstractVector{Symbol}}, i) = getindex(A, Not(I.skip...), i)
@inline Base.getindex(A::NamedRowMatrix, I::InvertedIndex{<:AbstractVector{Symbol}}, i::Integer) = getindex(A, Not(I.skip...), i)

@propagate_inbounds function Base.getindex(A::NamedRowVector, idx::InvertedIndex)
    idx_int = _inverted_symbols_to_ints(A, idx)
    return NamedRowVector(A.data[idx_int], reaxis(A, idx_int))
end
@propagate_inbounds function Base.getindex(A::NamedRowMatrix, idx, I::InvertedIndex)
    throw(ArgumentError("only named rows are allowed, not columns"))
end
@propagate_inbounds function Base.getindex(A::NamedRowMatrix, idx::InvertedIndex, I::InvertedIndex)
    throw(ArgumentError("only named rows are allowed, not columns"))
end
@propagate_inbounds function Base.getindex(A::NamedRowMatrix, I::InvertedIndex{<:AbstractVector{Symbol}}, i::InvertedIndex)
    throw(ArgumentError("only named rows are allowed, not columns"))
end
@propagate_inbounds function Base.getindex(A::NamedRowMatrix, idx::InvertedIndex, I)
    idx_int = _inverted_symbols_to_ints(A, idx)
    return NamedRowMatrix(A.data[idx_int, I], reaxis(A, idx_int))
end
@propagate_inbounds function Base.getindex(A::NamedRowMatrix, idx::InvertedIndex, I::Integer)
    idx_int = _inverted_symbols_to_ints(A, idx)
    return NamedRowVector(A.data[idx_int, I], reaxis(A, idx_int))
end
@propagate_inbounds function Base.getindex(A::NamedRowMatrix, idx::InvertedIndex{<:Union{Symbol,InvertedIndices.NotMultiIndex}})
    idx_int = _inverted_symbols_to_ints(A, idx)
    return NamedRowMatrix(A.data[idx_int, (:)], reaxis(A, idx_int))
end

# inverted indices
@inline Base.view(A::NamedRowVector, I::InvertedIndex{<:AbstractVector{Symbol}}) = view(A, Not(I.skip...))
@inline Base.view(A::NamedRowMatrix, I::InvertedIndex{<:AbstractVector{Symbol}}) = view(A, Not(I.skip...))
@inline Base.view(A::NamedRowMatrix, I::InvertedIndex{<:AbstractVector{Symbol}}, i) = view(A, Not(I.skip...), i)
@inline Base.view(A::NamedRowMatrix, I::InvertedIndex{<:AbstractVector{Symbol}}, i::Integer) = view(A, Not(I.skip...), i)

@propagate_inbounds function Base.view(A::NamedRowVector, idx::InvertedIndex)
    idx_int = _inverted_symbols_to_ints(A, idx)
    return NamedRowVector(view(A.data, idx_int), reaxis(A, idx_int))
end
@propagate_inbounds function Base.view(A::NamedRowMatrix, idx, I::InvertedIndex)
    throw(ArgumentError("only named rows are allowed, not columns"))
end
@propagate_inbounds function Base.view(A::NamedRowMatrix, idx::InvertedIndex, I::InvertedIndex)
    throw(ArgumentError("only named rows are allowed, not columns"))
end
@propagate_inbounds function Base.view(A::NamedRowMatrix, I::InvertedIndex{<:AbstractVector{Symbol}}, i::InvertedIndex)
    throw(ArgumentError("only named rows are allowed, not columns"))
end
@propagate_inbounds function Base.view(A::NamedRowMatrix, idx::InvertedIndex, I)
    idx_int = _inverted_symbols_to_ints(A, idx)
    return NamedRowMatrix(view(A.data, idx_int, I), reaxis(A, idx_int))
end
@propagate_inbounds function Base.view(A::NamedRowMatrix, idx::InvertedIndex, I::Integer)
    idx_int = _inverted_symbols_to_ints(A, idx)
    return NamedRowVector(view(A.data, idx_int, I), reaxis(A, idx_int))
end
@propagate_inbounds function Base.view(A::NamedRowMatrix, idx::InvertedIndex{<:Union{Symbol,InvertedIndices.NotMultiIndex}})
    idx_int = _inverted_symbols_to_ints(A, idx)
    return NamedRowMatrix(view(A.data, idx_int, (:)), reaxis(A, idx_int))
end

_inverted_symbols_to_ints(A, idx) = idx
function _inverted_symbols_to_ints(A, idx::InvertedIndex{Symbol})
    try
        Not(to_index(A, idx.skip))
    catch
        (:)
    end
end
function _inverted_symbols_to_ints(A, idx::InvertedIndex{InvertedIndices.NotMultiIndex}) 
    skips = Int64[]
    for i in idx.skip.indices
        try
            ind = to_index(A, i)
            push!(skips, ind)
        catch
        end
    end
    return Not(skips)
end


function Base.to_index(A::NamedRowArray, idx::Symbol)
    ind = findfirst(==(idx), A.rownames)
    isnothing(ind) && throw(ArgumentError("invalid index: $idx"))
    return ind
end

function Base.to_index(A::NamedRowArray, I::InvertedIndices.NotMultiIndex)
    [to_index(A, idx) for idx in I.indices]
end

function Base.to_indices(A::NamedRowArray, inds, I::Tuple{InvertedIndex{InvertedIndices.NotMultiIndex}, Vararg{Any}})
    new_indices = to_indices(A, inds, (I[1].skip, Base.tail(I)...))
    skips = InvertedIndices.uniquesort(new_indices[1])
    picks = InvertedIndices.spanned_indices(inds, skips)[1]
    return (InvertedIndices.InvertedIndexIterator(skips, picks), Base.tail(new_indices)...)
end

# Cartesian iteration
Base.eachindex(A::NamedRowArray) = eachindex(A.data)

"""
    Base.names(A::NamedRowArray)

Returns the row names of A as an NTuple.
"""
Base.names(A::NamedRowArray) = A.rownames


# Operator behavior

# unary operators
Base.:-(A::NamedRowArray) = NamedRowArray(-A.data, A.rownames)

# binary operators
Base.:*(a::Number, A::NamedRowArray) = NamedRowArray(a * A.data, A.rownames)
Base.:*(a::AbstractMatrix, A::NamedRowArray) = NamedRowArray(a * A.data, A.rownames)
Base.:*(A::NamedRowArray, a::Number) = NamedRowArray(A.data * a, A.rownames)
Base.:/(A::NamedRowArray, a::Number) = NamedRowArray(A.data / a, A.rownames)


# Broadcasting behavior (from https://docs.julialang.org/en/v1/manual/interfaces/#man-interfaces-broadcasting)
Base.BroadcastStyle(::Type{<:NamedRowArray}) = Broadcast.ArrayStyle{NamedRowArray}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{NamedRowArray}}, ::Type{ElType}) where ElType
    # Scan the inputs for the NamedRowArray:
    A = _find_nra(bc)
    # create the output
    NamedRowArray(similar(Array{ElType}, axes(bc)), A.rownames)
end

"`A = find_aac(As)` returns the first NamedRowArray among the arguments."
_find_nra(bc::Base.Broadcast.Broadcasted) = _find_nra(bc.args)
_find_nra(args::Tuple) = _find_nra(_find_nra(args[1]), Base.tail(args))
_find_nra(x) = x
_find_nra(::Tuple{}) = nothing
_find_nra(a::NamedRowArray, rest) = a
_find_nra(::Any, rest) = _find_nra(rest)

end # module
