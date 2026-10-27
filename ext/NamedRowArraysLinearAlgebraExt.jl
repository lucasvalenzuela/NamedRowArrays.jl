module NamedRowArraysLinearAlgebraExt

using NamedRowArrays: NamedRowArrays, NamedRowMatrix, NamedRowVector
using LinearAlgebra: Adjoint, Diagonal, Transpose

# The rows of a product are the rows of the left factor.
Base.:*(A::NamedRowMatrix, B::AbstractVector) = NamedRowArrays._withnames(A.data * B, A.rownames)
Base.:*(A::NamedRowMatrix, B::AbstractMatrix) = NamedRowArrays._withnames(A.data * B, A.rownames)
Base.:*(A::NamedRowMatrix, B::NamedRowVector) = NamedRowArrays._withnames(A.data * B.data, A.rownames)
Base.:*(A::NamedRowMatrix, B::NamedRowMatrix) = NamedRowArrays._withnames(A.data * B.data, A.rownames)

# Without names on the left factor, the product is a plain array. The data are passed on
# so that BLAS can be used for them.
Base.:*(A::AbstractMatrix, B::NamedRowVector) = A * B.data
Base.:*(A::AbstractMatrix, B::NamedRowMatrix) = A * B.data

# Methods resolving ambiguities with LinearAlgebra
Base.:*(A::Adjoint{<:Any,<:AbstractVector}, B::NamedRowMatrix) = A * B.data
Base.:*(A::Transpose{<:Any,<:AbstractVector}, B::NamedRowMatrix) = A * B.data
Base.:*(A::Union{Adjoint{T,V},Transpose{T,V}}, B::NamedRowVector) where {T,V<:AbstractVector} = A * B.data
Base.:*(A::Adjoint{<:Number,<:AbstractVector}, B::NamedRowVector{<:Number}) = A * B.data
Base.:*(A::Transpose{T,<:AbstractVector}, B::NamedRowVector{T}) where {T<:Real} = A * B.data
Base.:*(A::Diagonal, B::NamedRowVector) = A * B.data

end
