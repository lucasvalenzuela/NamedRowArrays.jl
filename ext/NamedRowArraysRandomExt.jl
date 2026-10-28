module NamedRowArraysRandomExt

using NamedRowArrays: NamedRowArrays, NamedRowMatrix, NamedRowVector, NamedRowArray
using Random: Random, AbstractRNG, randperm

Random.shuffle(rng::AbstractRNG, A::NamedRowVector) = A[NamedRowArrays._rowrange(A)[randperm(rng, length(A))]]
Random.shuffle(::AbstractRNG, ::NamedRowMatrix) = throw(
    ArgumentError(
        "shuffling a NamedRowMatrix moves elements between rows, so the rows would no longer match their names; " *
        "use A[shuffle(axes(A, 1)), :] to shuffle whole rows",
    ),
)

Random.shuffle!(::AbstractRNG, ::NamedRowArray) =
    NamedRowArrays._throw_inplace_reorder(Random.shuffle!, "shuffle or A[shuffle(axes(A, 1)), :]")

end
