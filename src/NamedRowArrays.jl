module NamedRowArrays

using Base: @propagate_inbounds, to_index
using InvertedIndices

export NamedRowArray, NamedRowVector, NamedRowMatrix, Not
public rownames

struct _Unchecked end

_lengtherror(nnames, nrows) =
    ArgumentError("the number of row names ($nnames) does not match the number of rows ($nrows)")

# e.g., ":a appears" or ":a, :c appear"
_appear(names) = length(names) == 1 ? "$(repr(only(names))) appears" : "$(join(repr.(names), ", ")) appear"

"""
    NamedRowArray(A::AbstractVecOrMat, names)

Wrap the vector or matrix `A` with a name for each row, so that rows can be selected by name as
well as by position. `A` is not copied. The row `names` are a collection of unique `Symbol`s, such
as a `Tuple` or a `Vector`, and are copied.

[`NamedRowVector{T,D,R}`](@ref NamedRowVector) and [`NamedRowMatrix{T,D,R}`](@ref NamedRowMatrix)
are aliases for `NamedRowArray{T,1,D,R}` and `NamedRowArray{T,2,D,R}`, where `T` is the element
type, `D` the type of the wrapped array and `R` the type of the vector of row names.

# Indexing

A `NamedRowArray` can be indexed like the wrapped array. In addition, the first index can be a row
name, a vector of row names, or a `Not` of row names, which selects all other rows. Names in
`Not` that the array does not have are ignored. A single row name, vector of row names or `Not`
of row names as the only index of a matrix selects whole rows. Row names in any other index
position throw an error.

The result is a `NamedRowArray` if all dimensions are indexed and the first index selects a
vector of rows; otherwise, it is the same as indexing the wrapped array. An index that selects
a row more than once throws an error, because the row names must be unique.

# Operations

- Broadcasting, `map`, arithmetic with numbers, reductions that keep every row (e.g.,
  `sum(A; dims=2)`) and `hcat` return a `NamedRowArray`. All `NamedRowArray`s among the
  arguments of broadcasting and `hcat` must have the same row names.
- The product `A * B` of a `NamedRowMatrix` `A` has the row names of `A`. If only `B` has row
  names, the product is a plain array.
- [`reverse`](@extref `Base.reverse-Tuple{AbstractVector}`), [`sort`](@extref Base.sort),
  [`sortslices`](@extref Base.sortslices), [`circshift`](@extref Base.circshift) and
  [`Random.shuffle`](@extref) reorder the row names together with the rows. Their in-place
  versions throw an error if they would reorder the rows, because other arrays can share the
  data (for example, views) and their row names would not be reordered. Reordering within the
  rows, e.g., `sort!(A; dims=2)`, is possible.
- `==` and `isequal` compare the values and ignore the row names.

# Examples

```jldoctest
julia> pos = NamedRowArray([1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0], (:x, :y, :z));

julia> pos[:y, 2]
5.0

julia> pos[:y]
3-element Vector{Float64}:
 4.0
 5.0
 6.0

julia> pos[2:3, 1]
2-element NamedRowVector(::Vector{Float64}) with eltype Float64 with row names:
 y z
 4.0
 7.0

julia> pos[[:z, :x], 2:3]
2×2 NamedRowMatrix(::Matrix{Float64}) with eltype Float64 with row names:
 z x
 8.0  9.0
 2.0  3.0

julia> pos[Not(:y), 2]
2-element NamedRowVector(::Vector{Float64}) with eltype Float64 with row names:
 x z
 2.0
 8.0
```
"""
struct NamedRowArray{T,N,D<:AbstractArray{T,N},R<:AbstractVector{Symbol}} <: AbstractArray{T,N}
    data::D
    # Views share the row names with their parent, so row names must never be modified
    # after construction.
    rownames::R

    function NamedRowArray{T,N,D}(data, names) where {T,N,D}
        N in (1, 2) ||
            throw(ArgumentError("a NamedRowArray wraps a vector or matrix, got a $N-dimensional array"))
        all(name -> name isa Symbol, names) || throw(ArgumentError("the row names must be Symbols"))
        # a copy, so that the row names cannot be changed from outside
        v = collect(Symbol, names)
        length(v) == size(data, 1) || throw(_lengtherror(length(v), size(data, 1)))
        allunique(v) || throw(
            ArgumentError("the row names must be unique, but $(_appear(_duplicates(v))) more than once"),
        )
        return new{T,N,D,Vector{Symbol}}(data, v)
    end

    # For row names that are known to be unique, such as those of another NamedRowArray.
    # They are stored as given.
    function NamedRowArray{T,N,D,R}(data, names, ::_Unchecked) where {T,N,D,R}
        length(names) == size(data, 1) || throw(_lengtherror(length(names), size(data, 1)))
        return new{T,N,D,R}(data, names)
    end
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

NamedRowArray(A::AbstractArray{T,N}, names) where {T,N} = NamedRowArray{T,N,typeof(A)}(A, names)
NamedRowVector(A::AbstractVector, names) = NamedRowArray(A, names)
NamedRowMatrix(A::AbstractMatrix, names) = NamedRowArray(A, names)

# With a copy of the row names of another NamedRowArray
_withnames(data, names) = _withownnames(data, collect(names))
# With row names that nothing else refers to
_withownnames(data::AbstractArray{T,N}, names::Vector{Symbol}) where {T,N} =
    NamedRowArray{T,N,typeof(data),Vector{Symbol}}(data, names, _Unchecked())
# With row names shared with another NamedRowArray
_withsharednames(data::AbstractArray{T,N}, names::AbstractVector{Symbol}) where {T,N} =
    NamedRowArray{T,N,typeof(data),typeof(names)}(data, names, _Unchecked())

_duplicates(names) = unique(filter(r -> count(==(r), names) > 1, names))

"""
    NamedRowArrays.rownames(A::NamedRowArray)

Return a copy of the row names of `A` as a `Vector{Symbol}`.
"""
rownames(A::NamedRowArray) = copy(A.rownames)

@inline Base.size(A::NamedRowArray) = size(A.data)
@inline Base.axes(A::NamedRowArray) = axes(A.data)
Base.IndexStyle(::Type{<:NamedRowArray{T,N,D}}) where {T,N,D} = IndexStyle(D)
Base.convert(::Type{Array{T,N}}, A::NamedRowArray{T,N}) where {T,N} = convert(Array{T,N}, A.data)
Base.parent(A::NamedRowArray) = A.data
Base.similar(A::NamedRowArray, ::Type{S}) where {S} = _withnames(similar(A.data, S), A.rownames)
Base.similar(A::NamedRowArray, ::Type{S}, dims::Dims{N}) where {S,N} = similar(A.data, S, dims)
function Base.reinterpret(::Type{T}, A::NamedRowArray) where {T}
    S = eltype(A)
    sizeof(T) == sizeof(S) || throw(
        ArgumentError(
            "cannot reinterpret a NamedRowArray with element type $S ($(sizeof(S)) bytes) as $T ($(sizeof(T)) bytes): " *
            "the number of rows would change, so the row names would no longer match",
        ),
    )
    return _withnames(reinterpret(T, A.data), A.rownames)
end


# Base's `summary` describes the array with this, e.g., as
# "3×4 NamedRowMatrix(view(::Matrix{Float64}, [1, 3, 4], :)) with eltype Float64".
function Base.showarg(io::IO, A::NamedRowArray, toplevel)
    print(io, A isa NamedRowVector ? "NamedRowVector(" : "NamedRowMatrix(")
    Base.showarg(io, A.data, false)
    print(io, ")")
    toplevel && print(io, " with eltype ", eltype(A))
    return nothing
end

# The row names go on their own line below the summary; Base prints the data.
function Base.show(io::IO, mime::MIME"text/plain", A::NamedRowArray)
    summary(io, A)
    println(io, " with row names:")
    _print_rownames(io, A.rownames)
    isempty(A) && return
    height, width = displaysize(io)
    data = sprint(show, mime, A.data; context=IOContext(io, :displaysize => (height - 1, width)))
    # drop the summary line of the data
    print(io, data[something(findfirst('\n', data), lastindex(data) + 1):end])
end

function _print_rownames(io, names::AbstractVector{Symbol})
    width = displaysize(io)[2]
    str = join(names, " ")
    if length(str) > width - 1
        # cut after the last complete name that fits, or inside a name wider than the display
        cut = first(str, max(width - 2, 0))
        ind = findlast(' ', cut)
        print(io, " ", isnothing(ind) ? cut : cut[1:ind], "…")
    else
        print(io, " ", str)
    end
end


# Indexing. Row names in the first index are converted to indices along the first axis of
# the data, and then all indices are converted by `to_indices` of the data. The result keeps
# the row names if the array is indexed in all of its dimensions and the row index selects a
# vector of rows.

@propagate_inbounds Base.getindex(A::NamedRowArray, I::Int...) = A.data[I...]
@propagate_inbounds function Base.getindex(A::NamedRowArray, I...)
    J = _to_indices(A, I)
    return _named(A, A.data[J...], J)
end

@propagate_inbounds function Base.view(A::NamedRowArray, I...)
    J = _to_indices(A, I)
    return _namedview(A, view(A.data, J...), I, J)
end

@propagate_inbounds function Base.setindex!(A::NamedRowArray, v, I::Int...)
    A.data[I...] = v
    return A
end
# The data converts the remaining indices itself: `setindex!` of an `Array` does not accept
# every index type that `to_indices` returns, e.g., `Base.LogicalIndex` on Julia 1.14.
@propagate_inbounds function Base.setindex!(A::NamedRowArray, v, I...)
    A.data[_rowindices(A, I)...] = v
    return A
end

_to_indices(A::NamedRowArray, I::Tuple) = to_indices(A.data, _rowindices(A, I))

# A single index of row names selects whole rows of a matrix.
_rowindices(A::NamedRowMatrix, I::Tuple{Any}) = _isrowname(I[1]) ? _rowindices(A, (I[1], :)) : I
function _rowindices(A::NamedRowArray, I::Tuple)
    any(_isrowname, Base.tail(I)) && throw(ArgumentError("only named rows are allowed, not columns"))
    return (_rowindex(A, I[1]), Base.tail(I)...)
end
_rowindices(A::NamedRowArray, I::Tuple{}) = ()

_selects_unique(_) = false
_selects_unique(::AbstractRange) = true
_selects_unique(::Base.LogicalIndex) = true
_selects_unique(::InvertedIndices.InvertedIndexIterator) = true

_isrowname(_) = false
_isrowname(::Symbol) = true
_isrowname(::AbstractVector{Symbol}) = true
_isrowname(I::InvertedIndex) = _isrowname(I.skip)
_isrowname(I::InvertedIndices.NotMultiIndex) = any(_isrowname, I.indices)

_keepsnames(A, J) = length(J) == ndims(A) && J[1] isa AbstractVector{<:Integer}

function _named(A, data, J)
    _keepsnames(A, J) || return data
    rows = J[1]
    rows isa Base.Slice && return _withnames(data, A.rownames)
    return _withownnames(data, _selectednames(A, rows))
end

function _namedview(A, data, I, J)
    _keepsnames(A, J) || return data
    rows = J[1]
    (first(axes(A.data, 1)) != 1 || _isintegervector(I[1])) &&
        return _withownnames(data, _selectednames(A, rows))
    # Reuse the row index of the data view, which Base creates for `Not` and logical indices.
    if data isa SubArray && parent(data) === A.data
        names = view(A.rownames, parentindices(data)[1])
    else
        names = view(A.rownames, rows)
    end
    _checkrepeats(rows, names)
    return _withsharednames(data, names)
end

# Integer vectors given by the caller can be changed after creating a view, unlike ranges,
# and unlike the index vectors that Base creates for logical indices.
_isintegervector(I) = I isa AbstractArray{<:Integer} && !(I isa AbstractRange) && !(I isa AbstractArray{Bool})

# A new vector of the row names that `rows`, indices along the first axis of the data, select
function _selectednames(A, rows)
    offset = first(axes(A.data, 1)) - 1
    if offset == 0
        names = A.rownames[rows]
    else
        names = Vector{Symbol}(undef, length(rows))
        for (k, i) in enumerate(rows)
            names[k] = A.rownames[i - offset]
        end
    end
    _checkrepeats(rows, names)
    return names
end

# only an index vector can select a row more than once
function _checkrepeats(rows, names)
    _selects_unique(rows) ||
        allunique(names) ||
        throw(
            ArgumentError(
                "the index selects $(join(repr.(_duplicates(names)), ", ")) more than once, " *
                "but the row names of a NamedRowArray must be unique",
            ),
        )
    return nothing
end

# Index along the first axis of the data of the row `name`, or `nothing` if there is no such row.
function _findrow(A::NamedRowArray, name::Symbol)
    pos = findfirst(==(name), A.rownames)
    return isnothing(pos) ? nothing : first(axes(A.data, 1)) + pos - 1
end

_rowindex(A, I) = I
function _rowindex(A, name::Symbol)
    ind = _findrow(A, name)
    isnothing(ind) && throw(ArgumentError("invalid row name: $name"))
    return ind
end
_rowindex(A, names::AbstractVector{Symbol}) = [_rowindex(A, name) for name in names]
_rowindex(A, I::InvertedIndex) = _isrowname(I.skip) ? Not(_skippedrows(A, I.skip)) : I
_rowindex(A, I::InvertedIndex{InvertedIndices.NotMultiIndex}) = Not(_skippedrows(A, I.skip))
function _rowindex(A, I::InvertedIndex{Symbol})
    ind = _findrow(A, I.skip)
    # the same type whether the name exists or not
    return isnothing(ind) ? Not(1:0) : Not(ind:ind)
end

# Row names that do not exist are ignored, so that `Not` can exclude names
# regardless of whether the array has them.
_skippedrows(A, I) = vec(collect(Int, to_index(A.data, I)))
function _skippedrows(A, name::Symbol)
    ind = _findrow(A, name)
    return isnothing(ind) ? Int[] : [ind]
end
_skippedrows(A, names::Union{AbstractVector{Symbol},InvertedIndices.NotMultiIndex}) =
    _skippedrows!(Int[], A, names)
function _skippedrows!(skips, A, names::AbstractVector{Symbol})
    for name in names
        ind = _findrow(A, name)
        isnothing(ind) || push!(skips, ind)
    end
    return skips
end
function _skippedrows!(skips, A, I::InvertedIndices.NotMultiIndex)
    for i in I.indices
        append!(skips, _skippedrows(A, i))
    end
    return skips
end

# Arithmetic with numbers keeps the type of the wrapped array, e.g., a lazy array, which
# broadcasting would turn into an Array.
Base.:-(A::NamedRowArray) = _withnames(-A.data, A.rownames)
Base.:*(a::Number, A::NamedRowArray) = _withnames(a * A.data, A.rownames)
Base.:*(A::NamedRowArray, a::Number) = _withnames(A.data * a, A.rownames)
Base.:/(A::NamedRowArray, a::Number) = _withnames(A.data / a, A.rownames)

_unnamed(A) = A
_unnamed(A::NamedRowArray) = A.data

# Mapping keeps the row names, as broadcasting does.
Base.map(f, A::NamedRowArray) = _withnames(map(f, A.data), A.rownames)

# Reductions that keep every row, e.g., `sum(A; dims=2)`, keep the row names.
function Base.mapreduce(f, op, A::NamedRowArray; dims=:, kwargs...)
    r = mapreduce(f, op, A.data; dims, kwargs...)
    _alongrows(dims) && return r
    return _withnames(r, A.rownames)
end

# Concatenating columns keeps the rows.
Base.hcat(A::NamedRowArray, B::AbstractVecOrMat...) = _hcat(A, B...)
Base.hcat(A::AbstractVecOrMat, B::NamedRowArray, C::AbstractVecOrMat...) = _hcat(A, B, C...)
Base.hcat(A::NamedRowArray, B::NamedRowArray, C::AbstractVecOrMat...) = _hcat(A, B, C...)
_hcat(As...) = _withnames(hcat(map(_unnamed, As)...), _commonrownames(As))


# Reordering. The data of a NamedRowArray can be shared with other arrays, e.g., views, whose
# row names would not follow if the rows were reordered in place. Functions that reorder rows
# in place therefore throw an error.

# Indices along the first axis of the data, as a range indexed from 1
_rowrange(A::NamedRowArray) = first(axes(A.data, 1)):last(axes(A.data, 1))

_alongrows(dims) = dims isa Colon || 1 in dims

function _throw_inplace_reorder(f, alternative)
    throw(
        ArgumentError(
            "$f cannot reorder the row names of a NamedRowArray in place; use $alternative instead, " *
            "which returns a new NamedRowArray with the row names reordered",
        ),
    )
end

function Base.reverse(A::NamedRowVector, start::Integer, stop::Integer=lastindex(A))
    r = _rowrange(A)
    return A[[first(r):(start - 1); stop:-1:start; (stop + 1):last(r)]]
end
Base.reverse(A::NamedRowVector; dims=:) = _reverse(A, dims)
Base.reverse(A::NamedRowMatrix; dims=:) = _reverse(A, dims)
function _reverse(A, dims)
    names = _alongrows(dims) ? reverse(A.rownames) : collect(A.rownames)
    return _withownnames(reverse(A.data; dims), names)
end

Base.reverse!(::NamedRowVector, ::Integer, ::Integer...) = _throw_inplace_reorder(reverse!, "reverse")
Base.reverse!(::NamedRowVector; dims=:) = _throw_inplace_reorder(reverse!, "reverse")
function Base.reverse!(A::NamedRowMatrix; dims=:)
    _alongrows(dims) && _throw_inplace_reorder(reverse!, "reverse")
    reverse!(A.data; dims)
    return A
end

Base.sort(A::NamedRowVector; kwargs...) = A[sortperm(A.data; kwargs...)]
function Base.sort(A::NamedRowMatrix; dims::Integer, kwargs...)
    dims == 1 && throw(
        ArgumentError(
            "sorting a NamedRowMatrix along dims=1 sorts each column separately, so the rows would no longer " *
            "match their names; use sortslices(A; dims=1) to reorder whole rows",
        ),
    )
    return _withnames(sort(A.data; dims, kwargs...), A.rownames)
end

Base.sort!(::NamedRowVector; kwargs...) = _throw_inplace_reorder(sort!, "sort")
function Base.sort!(A::NamedRowMatrix; dims::Integer, kwargs...)
    dims == 1 && _throw_inplace_reorder(sort!, "sortslices(A; dims=1)")
    sort!(A.data; dims, kwargs...)
    return A
end

function Base.sortslices(A::NamedRowMatrix; dims, kwargs...)
    dims == 2 && return _withnames(sortslices(A.data; dims, kwargs...), A.rownames)
    dims == 1 || throw(ArgumentError("dims must be 1 or 2 for a NamedRowMatrix, got $dims"))
    r = _rowrange(A)
    rows = [view(A.data, i, :) for i in r]
    return A[r[sortperm(rows; kwargs...)], :]
end

Base.permute!(::NamedRowVector, p::AbstractVector) = _throw_inplace_reorder(permute!, "A[p]")
Base.invpermute!(::NamedRowVector, p::AbstractVector) = _throw_inplace_reorder(invpermute!, "A[invperm(p)]")

Base.circshift(A::NamedRowArray, shiftamt::Real) = _circshift(A, shiftamt)
Base.circshift(A::NamedRowArray, shiftamt::Tuple{Vararg{Integer}}) = _circshift(A, shiftamt)
Base.circshift(A::NamedRowArray, shiftamt::AbstractVector{<:Integer}) = _circshift(A, shiftamt)
_circshift(A, shiftamt) =
    _withownnames(circshift(A.data, shiftamt), circshift(A.rownames, _rowshift(shiftamt)))

# A shift along the rows changes which values belong to the fixed row names of `dest`.
function _check_rowshift(dest, rowshift)
    mod(rowshift, size(dest, 1)) == 0 || _throw_inplace_reorder(circshift!, "circshift")
end
_rowshift(shiftamt::Real) = shiftamt
_rowshift(shiftamt) = isempty(shiftamt) ? 0 : first(shiftamt)

function Base.circshift!(dest::NamedRowArray, src, shiftamt)
    _check_rowshift(dest, _rowshift(shiftamt))
    circshift!(dest.data, _unnamed(src), shiftamt)
    return dest
end
Base.circshift!(dest::NamedRowArray, src, shiftamt::Tuple{}) =
    (circshift!(dest.data, _unnamed(src), shiftamt); dest)
function Base.circshift!(dest::NamedRowArray, src, shiftamt::Tuple{Vararg{Integer}})
    _check_rowshift(dest, _rowshift(shiftamt))
    circshift!(dest.data, _unnamed(src), shiftamt)
    return dest
end
function Base.circshift!(A::NamedRowVector, shift::Integer)
    _check_rowshift(A, shift)
    circshift!(A.data, shift)
    return A
end


# Broadcasting
Base.BroadcastStyle(::Type{<:NamedRowArray}) = Broadcast.ArrayStyle{NamedRowArray}()

function Base.similar(
    bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{NamedRowArray}},
    ::Type{ElType},
) where {ElType}
    _withnames(similar(Array{ElType}, axes(bc)), _commonrownames(bc))
end

# Row names of all NamedRowArrays among broadcast or concatenation arguments, which must agree.
_commonrownames(bc::Base.Broadcast.Broadcasted) = _commonrownames(bc.args)
_commonrownames(args::Tuple) =
    _combine_rownames(_commonrownames(first(args)), _commonrownames(Base.tail(args)))
_commonrownames(::Tuple{}) = nothing
_commonrownames(A::NamedRowArray) = A.rownames
_commonrownames(_) = nothing

_combine_rownames(::Nothing, ::Nothing) = nothing
_combine_rownames(a, ::Nothing) = a
_combine_rownames(::Nothing, b) = b
function _combine_rownames(a, b)
    a == b || throw(ArgumentError("cannot combine NamedRowArrays with different row names: $a and $b"))
    return a
end

end # module
