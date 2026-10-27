using NamedRowArrays
using Test
using Aqua
using ExplicitImports
using JET

using FillArrays
using LinearAlgebra
using OffsetArrays
using Random

fill4(v::AbstractVector, idx) = fill(4, count(!=(false), idx))
fill4(v::AbstractVector, idx::Colon) = fill(4, length(v))

fill4(a::AbstractMatrix, idx, idx2) = fill(4, count(!=(false), idx), count(!=(false), idx2))
fill4(a::AbstractMatrix, idx::Colon, idx2) = fill(4, size(a, 1), count(!=(false), idx2))
fill4(a::AbstractMatrix, idx, idx2::Colon) = fill(4, count(!=(false), idx), size(a, 2))
fill4(a::AbstractMatrix, idx::Colon, idx2::Colon) = fill(4, size(a, 1), size(a, 2))

@testset "NamedRowArrays.jl" begin
    @testset "Code quality" begin
        Aqua.test_all(NamedRowArrays)
        @test check_no_stale_explicit_imports(NamedRowArrays) === nothing
        @test check_no_self_qualified_accesses(NamedRowArrays) === nothing
        if JET.JET_AVAILABLE
            JET.test_package(NamedRowArrays; target_modules=(NamedRowArrays,))
        end
    end

    n = 4
    m = 5
    v = [
        0.6293451231426089
        0.4503389405961936
        0.47740714343281776
        0.7031298490032014
    ]
    a = [
        0.527015 0.928832 0.211228 0.933731 0.663758
        0.483022 0.830517 0.183429 0.365958 0.492058
        0.390663 0.666519 0.347508 0.56952 0.713226
        0.802763 0.133315 0.622647 0.0114776 0.165461
    ]
    rownames = Symbol.('a':'z')[1:n]

    arrs = [1:n, rand(m, n)', transpose(rand(m, n)), Fill(1.0, 4, m), @view rand(2n, m)[1:n, :]]

    idxs = [(:), 1:n, 2:3, 3:3, n:-1:1, [1, 3], [3, 1], [3], isodd.(1:n), collect(isodd.(1:n))]

    idxs2 = [(:), 1:m, 2:3, 3:3, m:-1:1, [1, 3], [3, 1], [3], isodd.(1:m), collect(isodd.(1:m))]
    nsymbols = [Not(:a), Not(:a, :c), Not([:b, :d]), Not(:d, :b), Not([:c, :a]), Not(:foo, :b), Not(:foo)]
    nints = [Not(1), Not(1, 3), Not([2, 4]), Not(4, 2), Not([3, 1]), Not(2), (:)]

    @testset "Construction" begin
        vn = NamedRowArray(v, rownames)
        @test vn.data === v
        @test vn.data !== vn
        @test vn.data == vn
        @test all(vn.rownames .=== rownames)
        @test NamedRowVector(v, rownames) == vn
        @test NamedRowArray(v, Tuple(rownames)) == vn
        @test NamedRowVector(v, Tuple(rownames)) == vn

        an = NamedRowArray(a, rownames)
        @test an.data === a
        @test an.data !== an
        @test an.data == an
        @test all(an.rownames .=== rownames)
        @test NamedRowMatrix(a, rownames) == an
        @test NamedRowArray(a, Tuple(rownames)) == an
        @test NamedRowMatrix(a, Tuple(rownames)) == an

        msg = "the number of row names (3) does not match the number of rows (4)"
        @test_throws msg NamedRowArray(v, rownames[1:3])
        @test_throws msg NamedRowArray(a, rownames[1:3])
        @test_throws "the row names must be Symbols" NamedRowArray(v, String.(rownames))
        @test_throws "the row names must be Symbols" NamedRowArray(a, String.(rownames))
        msg = "a NamedRowArray wraps a vector or matrix, got a 3-dimensional array"
        @test_throws msg NamedRowArray(rand(4, 4, 4), rownames)
        @test NamedRowArray(Float64[], Symbol[]).rownames == Symbol[]
        @test NamedRowArray(v, rownames)[falses(n)].rownames == Symbol[]
        @test NamedRowArray(a, rownames)[falses(n), :].rownames == Symbol[]
        msg = "the row names must be unique, but :a appears more than once"
        @test_throws msg NamedRowArray(v, [:a, :b, :a, :c])
        @test_throws "but :a, :c appear more than once" NamedRowArray(a, (:a, :c, :a, :c))

        @test_throws MethodError NamedRowArray(1, rownames)
        @test_throws MethodError NamedRowVector(a, rownames)
        @test_throws MethodError NamedRowMatrix(v, rownames)

        @test NamedRowArrays.rownames(vn) == vn.rownames
        @test NamedRowArrays.rownames(an) == an.rownames

        # the row names can be given as any collection of Symbols and are stored as a copy
        @test NamedRowArray(v, Tuple(rownames)).rownames == rownames
        @test NamedRowArray(v, (s for s in rownames)).rownames == rownames
        @test NamedRowArray(v, Tuple(rownames)) isa NamedRowArray{Float64,1,Vector{Float64}}
        @test NamedRowVector(v, Tuple(rownames)) isa NamedRowVector{Float64,Vector{Float64}}
        @test NamedRowMatrix(a, Tuple(rownames)) isa NamedRowMatrix{Float64,Matrix{Float64}}
        names_in = copy(rownames)
        vn2 = NamedRowArray(v, names_in)
        names_in[1] = :changed
        @test vn2.rownames == rownames
        NamedRowArrays.rownames(vn2)[1] = :changed
        @test vn2.rownames == rownames
        @test_throws "the row names must be Symbols" NamedRowArray(v, "abcd")

        for arr in arrs
            @test_nowarn NamedRowArray(arr, rownames)
        end
    end

    @testset "Array interface" begin
        vn = NamedRowArray(v, rownames)
        an = NamedRowArray(a, rownames)

        @test length(vn) == length(v)
        @test size(vn) == size(v)
        @test size(vn, 1) == size(v, 1)
        @test axes(vn) == axes(v)
        @test axes(vn, 1) == axes(v, 1)
        @test eachindex(vn) == eachindex(v)

        @test parent(vn) === v
        @test convert(Array, vn) == v
        @test convert(Array, vn) isa Vector
        @test collect(vn) isa Vector
        @test similar(vn) isa NamedRowVector
        @test similar(vn) !== vn
        @test size(similar(vn)) == size(vn)
        @test eltype(similar(vn)) == eltype(vn)
        @test similar(vn).rownames == vn.rownames
        @test similar(vn, Int) isa NamedRowVector{Int}
        @test similar(vn, Int).rownames == vn.rownames
        @test reinterpret(Int64, vn) isa NamedRowVector{Int64}
        @test reinterpret(Int64, vn).rownames == vn.rownames
        @test reinterpret(Float64, reinterpret(Int64, vn)) == vn
        msg = "cannot reinterpret a NamedRowArray with element type Float64 (8 bytes) as Float32 (4 bytes)"
        @test_throws msg reinterpret(Float32, vn)

        @test copy(vn) !== vn
        @test copy(vn).data !== vn.data
        @test copy(vn).rownames == vn.rownames
        @test deepcopy(vn) !== vn


        @test length(an) == length(a)
        @test size(an) == size(a)
        @test size(an, 1) == size(a, 1)
        @test axes(an) == axes(a)
        @test axes(an, 1) == axes(a, 1)
        @test eachindex(an) == eachindex(a)

        @test parent(an) === a
        @test convert(Array, an) == a
        @test convert(Array, an) isa Matrix
        @test collect(an) isa Matrix
        @test similar(an) isa NamedRowMatrix
        @test similar(an) !== an
        @test size(similar(an)) == size(an)
        @test eltype(similar(an)) == eltype(an)
        @test similar(an).rownames == an.rownames
        @test similar(an, Int) isa NamedRowMatrix{Int}
        @test similar(an, Int).rownames == an.rownames
        @test reinterpret(Int64, an) isa NamedRowMatrix{Int64}
        @test_throws "the number of rows would change" reinterpret(Int32, an)

        @test IndexStyle(an) == IndexLinear()
        @test IndexStyle(NamedRowArray(view(a, [1, 2, 3, 4], 1:3), rownames)) == IndexCartesian()
        @test IndexStyle(NamedRowArray(a', Symbol.('a':'e'))) == IndexCartesian()

        @test copy(an) !== an
        @test copy(an).data !== an.data
        @test copy(an).rownames == an.rownames
        @test deepcopy(an) !== an
    end

    @testset "getindex and view" begin
        vn = NamedRowArray(v, rownames)
        an = NamedRowArray(a, rownames)

        @testset "Vector" begin
            # single elements
            @test vn[1] == v[1]
            @test vn[2] == v[2]
            @test_throws BoundsError vn[0]
            @test_throws BoundsError vn[n + 1]

            @test vn[rownames[1]] == v[1]
            @test vn[rownames[2]] == v[2]
            @test_throws "invalid row name: foo" vn[:foo]
            @test_throws ArgumentError vn["foo"]

            @test @view(vn[rownames[1]]) == @view(v[1])
            @test @view(vn[rownames[2]]) == @view(v[2])
            @test_throws "invalid row name: foo" @view(vn[:foo])
            @test_throws ArgumentError @view(vn["foo"])

            # multiple elements
            for idx in idxs
                @test vn[idx] == v[idx]
                @test vn[idx] isa NamedRowVector
                @test all(vn[idx].rownames .== rownames[idx])

                @test @view(vn[idx]) == v[idx]
                @test @view(vn[idx]) isa NamedRowVector
                @test @view(vn[idx]).data isa SubArray
                @test all(@view(vn[idx]).rownames .== rownames[idx])
            end

            @test vn[rownames[2:3]] == v[2:3]
            @test vn[rownames[2:3]] isa NamedRowVector
            @test vn[rownames[2:3]].rownames == vn.rownames[2:3]
            @test_throws "invalid row name: foo" vn[[:foo, :a]]
            @test_throws ArgumentError vn[["foo", "a"]]

            @test @view(vn[rownames[2:3]]) == v[2:3]
            @test @view(vn[rownames[2:3]]) isa NamedRowVector
            @test @view(vn[rownames[2:3]]).data isa SubArray
            @test @view(vn[rownames[2:3]]).rownames == vn.rownames[2:3]
            @test_throws "invalid row name: foo" @view(vn[[:foo, :a]])
            @test_throws ArgumentError @view(vn[["foo", "a"]])

            # inverted
            for (nsymbol, nint) in zip(nsymbols, nints)
                @test vn[nint] == v[nint]
                @test vn[nint] isa NamedRowVector
                @test collect(vn[nint].rownames) == collect(vn.rownames)[nint]

                @test vn[nsymbol] == v[nint]
                @test vn[nsymbol] isa NamedRowVector
                @test collect(vn[nsymbol].rownames) == collect(vn.rownames)[nint]

                @test @view(vn[nint]) == v[nint]
                @test @view(vn[nint]) isa NamedRowVector
                @test @view(vn[nint]).data isa SubArray
                @test collect(@view(vn[nint]).rownames) == collect(vn.rownames)[nint]

                @test @view(vn[nsymbol]) == v[nint]
                @test @view(vn[nsymbol]) isa NamedRowVector
                @test @view(vn[nsymbol]).data isa SubArray
                @test collect(@view(vn[nsymbol]).rownames) == collect(vn.rownames)[nint]
            end
        end

        @testset "Matrix" begin
            # single elements
            @test an[1] == a[1]
            @test an[2] == a[2]
            @test an[n + 1] == a[n + 1]
            @test_throws BoundsError an[0]
            @test_throws BoundsError an[length(an) + 1]

            @test an[1, 1] == a[1, 1]
            @test an[2, 3] == a[2, 3]
            @test an[n, 3] == a[n, 3]
            @test_throws BoundsError an[n + 1, 1]
            @test_throws BoundsError an[0, 1]

            @test an[CartesianIndex(1, 2)] == a[CartesianIndex(1, 2)]
            @test an[[CartesianIndex(1, 2), CartesianIndex(2, 2)]] ==
                  a[[CartesianIndex(1, 2), CartesianIndex(2, 2)]]

            @test an[rownames[1], 1] == a[1, 1]
            @test an[rownames[2], 3] == a[2, 3]
            @test_throws "invalid row name: foo" an[:foo, 1]
            @test_throws ArgumentError an["foo", 1]
            @test_throws "only named rows are allowed, not columns" an[rownames[1], rownames[1]]

            @test @view(an[rownames[1], 1]) == @view(a[1, 1])
            @test @view(an[rownames[2], 3]) == @view(a[2, 3])
            @test_throws "invalid row name: foo" @view(an[:foo, 1])
            @test_throws ArgumentError @view(an["foo", 1])
            @test_throws "only named rows are allowed, not columns" @view(an[rownames[1], rownames[1]])

            # multiple elements
            for idx in idxs
                for idx2 in idxs2
                    @test an[idx, idx2] == a[idx, idx2]
                    @test an[idx, idx2] isa NamedRowMatrix
                    @test all(an[idx, idx2].rownames .== rownames[idx])
                end

                @test an[idx, 2] == a[idx, 2]
                @test an[idx, 2] isa NamedRowVector
                @test all(an[idx, 2].rownames .== rownames[idx])
            end

            for idx2 in idxs2
                @test an[rownames[1], idx2] == a[1, idx2]
                @test an[rownames[2], idx2] == a[2, idx2]
                @test_throws "invalid row name: foo" an[:foo, idx2]
                @test_throws ArgumentError an["foo", idx2]

                @test an[rownames[2:3], idx2] == a[2:3, idx2]
                @test an[rownames[2:3], idx2] isa NamedRowMatrix
                @test an[rownames[2:3], idx2].rownames == an.rownames[2:3]
                @test_throws "invalid row name: foo" an[[:foo, :a], idx2]
                @test_throws ArgumentError an[["foo", "a"], idx2]

                @test @view(an[rownames[1], idx2]) == a[1, idx2]
                @test @view(an[rownames[2], idx2]) == a[2, idx2]
                @test_throws "invalid row name: foo" @view(an[:foo, idx2])
                @test_throws ArgumentError @view(an["foo", idx2])

                @test @view(an[rownames[2:3], idx2]) == a[2:3, idx2]
                @test @view(an[rownames[2:3], idx2]) isa NamedRowMatrix
                @test @view(an[rownames[2:3], idx2]).data isa SubArray
                @test @view(an[rownames[2:3], idx2]).rownames == an.rownames[2:3]
                @test_throws "invalid row name: foo" @view(an[[:foo, :a], idx2])
                @test_throws ArgumentError @view(an[["foo", "a"], idx2])
            end

            @test an[rownames[1]] == a[1, :]
            @test an[rownames[2]] == a[2, :]
            @test an[rownames[1]] isa Vector
            @test_throws "invalid row name: foo" an[:foo]
            @test_throws ArgumentError an["foo"]

            @test an[rownames[2:3], 1] == a[2:3, 1]
            @test an[rownames[2:3], 3] == a[2:3, 3]
            @test an[rownames[2:3], 1] isa NamedRowVector
            @test an[rownames[2:3], 1].rownames == an.rownames[2:3]
            @test_throws "invalid row name: foo" an[[:foo, :a], 1]
            @test_throws ArgumentError an[["foo", "a"], 1]

            @test @view(an[rownames[1]]) == a[1, :]
            @test @view(an[rownames[2]]) == a[2, :]
            @test @view(an[rownames[1]]) isa SubArray
            @test_throws "invalid row name: foo" @view(an[:foo])
            @test_throws ArgumentError @view(an["foo"])

            @test @view(an[rownames[2:3], 1]) == a[2:3, 1]
            @test @view(an[rownames[2:3], 3]) == a[2:3, 3]
            @test @view(an[rownames[2:3], 1]) isa NamedRowVector
            @test @view(an[rownames[2:3], 1]).data isa SubArray
            @test @view(an[rownames[2:3], 1]).rownames == an.rownames[2:3]
            @test_throws "invalid row name: foo" @view(an[[:foo, :a], 1])
            @test_throws ArgumentError @view(an[["foo", "a"], 1])

            # inverted
            for (nsymbol, nint) in zip(nsymbols, nints)
                for idx2 in idxs2
                    @test an[nint, idx2] == a[nint, idx2]
                    @test an[nint, idx2] isa NamedRowMatrix
                    @test collect(an[nint, idx2].rownames) == collect(an.rownames)[nint]

                    @test an[nsymbol, idx2] == a[nint, idx2]
                    @test an[nsymbol, idx2] isa NamedRowMatrix
                    @test collect(an[nsymbol, idx2].rownames) == collect(an.rownames)[nint]

                    @test @view(an[nint, idx2]) == a[nint, idx2]
                    @test @view(an[nint, idx2]) isa NamedRowMatrix
                    @test @view(an[nint, idx2]).data isa SubArray
                    @test collect(@view(an[nint, idx2]).rownames) == collect(an.rownames)[nint]

                    @test @view(an[nsymbol, idx2]) == a[nint, idx2]
                    @test @view(an[nsymbol, idx2]) isa NamedRowMatrix
                    @test @view(an[nsymbol, idx2]).data isa SubArray
                    @test collect(@view(an[nsymbol, idx2]).rownames) == collect(an.rownames)[nint]
                end

                @test an[nint] == a[nint]
                @test !(an[nint] isa NamedRowArray)

                @test an[nsymbol] == a[nint, :]
                @test an[nsymbol] isa NamedRowMatrix
                @test collect(an[nsymbol].rownames) == collect(an.rownames)[nint]

                @test an[nint, 2] == a[nint, 2]
                @test an[nint, 2] isa NamedRowVector

                @test an[nsymbol, 2] == a[nint, 2]
                @test an[nsymbol, 2] isa NamedRowVector
                @test collect(an[nsymbol, 2].rownames) == collect(an.rownames)[nint]

                @test_throws ArgumentError an[nsymbol, nsymbol]
                @test_throws ArgumentError an[1, nsymbol]
            end
        end
    end

    @testset "Index types" begin
        vn = NamedRowArray(v, rownames)
        an = NamedRowArray(a, rownames)

        # all integer types select rows by position
        @test vn[UInt[1, 3]] == v[[1, 3]]
        @test vn[UInt[1, 3]].rownames == [:a, :c]
        @test vn[Int32(2):Int32(3)].rownames == [:b, :c]
        @test an[UInt(2):UInt(3), 1].rownames == [:b, :c]
        @test view(an, Int8[4, 1], :).rownames == [:d, :a]
        @test vn[Int[]].rownames == Symbol[]

        # Cartesian and linear indices give plain values or arrays
        @test vn[CartesianIndex(2)] == v[2]
        @test an[CartesianIndex(2, 1)] == a[2, 1]
        @test an[[CartesianIndex(1, 1), CartesianIndex(2, 2)]] isa Vector
        @test an[2:3] isa Vector
        @test an[Not(1)] isa Vector
        @test vn[[1 2; 3 4]] == v[[1 2; 3 4]]
        @test vn[[1 2; 3 4]] isa Matrix

        # a single row index on a matrix gives plain arrays
        @test an[2, 1:2] isa Vector
        @test an[:b, 1:2] isa Vector

        @test_throws "invalid row name: foo" vn[:foo]
        @test_throws "invalid row name: foo" an[[:a, :foo], 1]
        msg = "only named rows are allowed, not columns"
        @test_throws msg an[1, :a]
        @test_throws msg an[:, [:a]]
        @test_throws msg an[:a, Not(:b)]
        @test_throws msg view(an, :, :a)
        @test_throws msg (copy(an)[1, :a] = 0)

        vnc = copy(vn)
        @test setindex!(vnc, 0, :b) === vnc
        @test setindex!(vnc, 0, 1) === vnc
    end

    @testset "Inverted indices" begin
        vn = NamedRowArray(v, rownames)
        # row names that do not exist are ignored
        @test vn[Not(:a, :foo)].rownames == [:b, :c, :d]
        @test vn[Not(:foo)].rownames == rownames
        @test vn[Not(:a, 3)].rownames == [:b, :d]
        # other invalid indices throw
        @test_throws ArgumentError vn[Not(:a, "b")]
    end

    @testset "Row names of views" begin
        an = NamedRowArray(a, rownames)
        shares(B, A) = B.rownames isa SubArray && parent(B.rownames) === A.rownames

        # views share the row names with their parent
        @test shares(view(an, 2:3, :), an)
        @test shares(view(an, :, 2), an)
        @test shares(view(an, Not(:a), 1), an)
        @test shares(view(an, Not(:a, :foo), 1:2), an)
        @test shares(view(an, [:c, :a], 1), an)
        @test shares(view(an, isodd.(1:n), :), an)
        @test shares(@view(an[[:b, :d]]), an)
        @test view(an, Not(:a), 1).rownames == [:b, :c, :d]
        @test view(an, [:c, :a], 1).rownames == [:c, :a]

        # views of views share the row names of the first array
        vv = view(view(an, 2:4, :), [:d, :b], 1)
        @test vv.rownames == [:d, :b]
        @test parent(vv.rownames) === an.rownames

        # an integer vector given by the caller can change, so the view gets its own row names
        idx = [3, 1]
        vi = view(an, idx, :)
        @test vi.rownames isa Vector{Symbol}
        idx[2] = 3
        @test vi.rownames == [:c, :a]

        # everything else has its own row names
        @test an[2:3, :].rownames isa Vector{Symbol}
        @test similar(view(an, 2:3, :)).rownames isa Vector{Symbol}
        @test NamedRowArrays.rownames(view(an, 2:3, :)) isa Vector{Symbol}
        @test (view(an, 2:3, :) .+ 1).rownames isa Vector{Symbol}
        @test view(NamedRowArray(OffsetArray(a, -1, 0), rownames), 1:2, :).rownames == [:b, :c]
        @test_throws "the index selects :a more than once" view(an, [:a, :b, :a], 1)
    end

    @testset "Repeated rows" begin
        vn = NamedRowArray(v, rownames)
        an = NamedRowArray(a, rownames)
        msg = "the index selects :a more than once, but the row names of a NamedRowArray must be unique"
        @test_throws msg vn[[1, 1, 2]]
        @test_throws msg vn[[:a, :b, :a]]
        @test_throws msg view(vn, [:a, :a])
        @test_throws "the index selects :a, :c more than once" an[[1, 3, 1, 3], :]
    end

    @testset "Inferred return types" begin
        vn = NamedRowArray(v, rownames)
        an = NamedRowArray(a, rownames)
        for I in [
            (2:3,),
            ([1, 3],),
            (:a,),
            ([:a, :c],),
            (Not(:a),),
            (Not(:a, :b),),
            (Not(:a, [2]),),
            (Not(:foo),),
            (Not(1),),
            (isodd.(1:n),),
            (:,),
        ]
            @inferred vn[I...]
            @inferred view(vn, I...)
        end
        for I in [
            (:a,),
            ([:a, :b],),
            (Not(:a),),
            (:, 2),
            (2:3, 1),
            (:a, 2:3),
            (2:3,),
            (Not(:b), :),
            (Not(:b, :foo), 1:2),
        ]
            @inferred an[I...]
            @inferred view(an, I...)
        end
    end

    @testset "setindex!" begin
        vn = NamedRowArray(v, rownames)
        an = NamedRowArray(a, rownames)

        @testset "Vector" begin
            # single elements
            vnc = copy(vn)
            vnc[1] = 4
            @test vnc[1] == 4
            vnc[2] = 4
            @test vnc[2] == 4
            @test_throws BoundsError (vnc[0] = 4)
            @test_throws BoundsError (vnc[n + 1] = 4)

            vnc = copy(vn)
            vnc[rownames[1]] = 4
            @test vnc[rownames[1]] == 4
            vnc[rownames[2]] = 4
            @test vnc[rownames[2]] == 4
            @test_throws "invalid row name: foo" (vnc[:foo] = 4)
            @test_throws ArgumentError (vnc["foo"] = 4)

            # multiple elements
            for idx in idxs
                vnc = copy(vn)
                vnc[idx] = fill4(vnc, idx)
                @test all(vnc[idx] .== 4)

                vnc = copy(vn)
                vnc[rownames[2:3]] = fill4(vnc, 2:3)
                @test all(vnc[rownames[2:3]] .== 4)
                @test_throws "invalid row name: foo" (vnc[[:foo, :a]] = fill4(vnc, 2:3))
                @test_throws ArgumentError (vnc[["foo", "a"]] = fill4(vnc, 2:3))
            end
        end

        @testset "Matrix" begin
            # single elements
            anc = copy(an)
            anc[1] = 4
            @test anc[1] == 4
            anc[2] = 4
            @test anc[2] == 4
            anc[n + 1] = 4
            @test anc[n + 1] == 4
            @test_throws BoundsError (anc[0] = 4)
            @test_throws BoundsError (anc[length(an) + 1] = 4)

            anc = copy(an)
            anc[1, 1] = 4
            @test anc[1, 1] == 4
            anc[2, 3] = 4
            @test anc[2, 3] == 4
            anc[n, 3] = 4
            @test anc[n, 3] == 4
            anc[CartesianIndex(3, 2)] = 4
            @test anc[CartesianIndex(3, 2)] == 4
            @test_throws BoundsError (anc[n + 1, 1] = 4)
            @test_throws BoundsError (anc[0, 1] = 4)

            anc = copy(an)
            anc[rownames[1], 1] = 4
            @test anc[rownames[1], 1] == 4
            anc[rownames[2], 3] = 4
            @test anc[rownames[2], 3] == 4
            @test_throws "invalid row name: foo" (anc[:foo, 1] = 4)
            @test_throws ArgumentError (anc["foo", 1] = 4)
            @test_throws "only named rows are allowed, not columns" (anc[rownames[1], rownames[1]] = 4)

            # multiple elements
            for idx in idxs
                for idx2 in idxs2
                    anc = copy(an)
                    anc[idx, idx2] = fill4(anc, idx, idx2)
                    @test all(anc[idx, idx2] .== 4)
                end

                anc = copy(an)
                anc[idx, 2] = fill4(anc, idx, 1)
                @test all(anc[idx, 2] .== 4)
            end

            for idx2 in idxs2
                anc = copy(an)
                anc[rownames[1], idx2] = fill4(anc, 1, idx2)
                @test all(anc[rownames[1], idx2] .== 4)
                anc[rownames[2], idx2] = fill4(anc, 1, idx2)
                @test all(anc[rownames[2], idx2] .== 4)
                @test_throws "invalid row name: foo" (anc[:foo, idx2] = fill4(anc, 1, idx2))
                @test_throws ArgumentError (anc["foo", idx2] = fill4(anc, 1, idx2))

                anc = copy(an)
                anc[rownames[2:3], idx2] = fill4(anc, 2:3, idx2)
                @test all(anc[rownames[2:3], idx2] .== 4)
                @test_throws "invalid row name: foo" (anc[[:foo, :a], idx2] = fill4(anc, 2:3, idx2))
                @test_throws ArgumentError (anc[["foo", "a"], idx2] = fill4(anc, 2:3, idx2))
            end

            anc = copy(an)
            anc[[CartesianIndex(3, 2), CartesianIndex(3, 3)]] = [4, 4]
            @test all(anc[[CartesianIndex(3, 2), CartesianIndex(3, 3)]] .== 4)

            anc = copy(an)
            anc[rownames[1]] = fill4(anc, 1, :)
            @test all(anc[rownames[1]] .== 4)
            anc[rownames[2]] = fill4(anc, 1, :)
            @test all(anc[rownames[2]] .== 4)
            @test_throws "invalid row name: foo" (anc[:foo] = 4)
            @test_throws ArgumentError (anc["foo"] = 4)

            anc = copy(an)
            anc[rownames[2:3]] = fill4(anc, 2:3, :)
            @test all(anc[rownames[2:3]] .== 4)
            @test_throws "invalid row name: foo" (anc[[:foo, :a]] = fill4(anc, 2:3, :))
            @test_throws ArgumentError (anc[["foo", "a"]] = fill4(anc, 2:3, :))
        end
    end

    @testset "Arrays with other axes" begin
        ov = NamedRowArray(OffsetArray([10, 20, 30, 40], 0:3), rownames)
        data = reshape(1:20, 4, 5)
        oa = NamedRowArray(OffsetArray(data, -1, 2), rownames)
        @test axes(ov) == (0:3,)

        @test ov[:a] == 10
        @test ov[:d] == 40
        @test ov[[:c, :a]] == [30, 10]
        @test ov[[:c, :a]].rownames == [:c, :a]
        @test ov[1:2].rownames == [:b, :c]
        @test ov[[3, 0]].rownames == [:d, :a]
        @test ov[OffsetArray([true, false, true, false], 0:3)].rownames == [:a, :c]
        @test ov[:].rownames == rownames
        @test view(ov, [:d, :b]) == [40, 20]
        @test view(ov, [:d, :b]).rownames == [:d, :b]
        @test ov[Not(:b)] == [10, 30, 40]
        @test ov[Not(:b)].rownames == [:a, :c, :d]
        @test ov[Not(:b, :foo, 3)].rownames == [:a, :c]
        @test ov[Not(0)].rownames == [:b, :c, :d]

        @test oa[:a, 3] == data[1, 1]
        @test parent(oa[:b]) == data[2, :]
        @test axes(oa[:b]) == (3:7,)
        @test oa[[:d, :a], 3:4] == data[[4, 1], 1:2]
        @test oa[[:d, :a], 3:4].rownames == [:d, :a]
        @test oa[1:2, 3].rownames == [:b, :c]
        @test oa[Not(:a), 4].rownames == [:b, :c, :d]

        ovc = copy(ov)
        ovc[:b] = 0
        @test ovc[1] == 0
        oac = copy(oa)
        oac[:c, 4] = 0
        @test oac[2, 4] == 0

        @test sort(ov; rev=true).rownames == [:d, :c, :b, :a]
        @test reverse(ov).rownames == [:d, :c, :b, :a]
        @test reverse(ov, 1, 2).rownames == [:a, :c, :b, :d]
        @test sortslices(oa; dims=1, rev=true).rownames == [:d, :c, :b, :a]
    end

    @testset "Broadcasting and arithmetic" begin
        vn = NamedRowArray(v, rownames)
        an = NamedRowArray(a, rownames)

        for arr in [vn, an]
            # operators
            @test +arr isa NamedRowArray
            @test -arr isa NamedRowArray
            @test 2arr isa NamedRowArray
            @test arr * 2 isa NamedRowArray
            @test arr / 2 isa NamedRowArray

            # broadcasting
            @test arr ./ 2 isa NamedRowArray
            @test arr .* 2 isa NamedRowArray
            @test exp10.(arr) .- 2 isa NamedRowArray
            @test log10.(arr) .+ 2 isa NamedRowArray
            @test log10.(arr) .+ arr isa NamedRowArray
            @test (arr .+ 2arr).rownames == arr.rownames
        end

        @test -vn == -v
        @test (-vn).rownames == vn.rownames
        @test 2an == 2a
        @test (an * 2).rownames == an.rownames
        @test (an / 2).rownames == an.rownames

        # arithmetic with numbers keeps the type of the wrapped array
        fn = NamedRowArray(Fill(2.0, n, m), rownames)
        for B in (-fn, 3fn, fn * 3, fn / 2)
            @test B isa NamedRowMatrix
            @test B.data isa Fill
            @test B.rownames == rownames
        end
        @test (fn / 2) == fill(1.0, n, m)

        @test (vn .+ an).rownames == an.rownames
        @test (an .+ v).rownames == an.rownames
        msg = "the number of row names (1) does not match the number of rows (4)"
        @test_throws msg NamedRowArray([1.0], [:a]) .+ v
        other = NamedRowArray(v, [:w, :x, :y, :z])
        msg = "cannot combine NamedRowArrays with different row names: [:a, :b, :c, :d] and [:w, :x, :y, :z]"
        @test_throws msg vn .+ other
        @test_throws "different row names" an .* 2 .+ exp.(other)

        # map keeps the row names and the type of the wrapped array
        @test map(x -> 2x, vn) == 2v
        @test map(x -> 2x, vn).rownames == rownames
        @test map(x -> 2x, an).rownames == rownames
        @test map(x -> 2x, NamedRowArray(Fill(1.0, n), rownames)).data isa Fill
    end

    @testset "Reductions" begin
        an = NamedRowArray(a, rownames)
        # reductions that keep every row keep the row names
        for f in (sum, prod, maximum, minimum)
            @test f(an; dims=2) == f(a; dims=2)
            @test f(an; dims=2).rownames == rownames
            @test f(an; dims=1) == f(a; dims=1)
            @test !(f(an; dims=1) isa NamedRowArray)
            @test f(an) == f(a)
        end
        @test sum(abs2, an; dims=2).rownames == rownames
        @test count(>(0.5), an; dims=2) == count(>(0.5), a; dims=2)
        @test count(>(0.5), an; dims=2).rownames == rownames
        @test any(>(0.5), an; dims=2).rownames == rownames
        @test !(sum(an; dims=(1, 2)) isa NamedRowArray)
        @test sum(NamedRowArray(Fill(1.0, n, m), rownames); dims=2) == fill(m, n, 1)
    end

    @testset "hcat" begin
        vn = NamedRowArray(v, rownames)
        an = NamedRowArray(a, rownames)
        for (B, expected) in (
            (hcat(an, an), hcat(a, a)),
            (hcat(vn, an), hcat(v, a)),
            (hcat(v, an), hcat(v, a)),
            (hcat(an, a, vn), hcat(a, a, v)),
            ([an v], [a v]),
            ([v vn an], [v v a]),
        )
            @test B == expected
            @test B isa NamedRowMatrix
            @test B.rownames == rownames
        end
        other = NamedRowArray(v, [:w, :x, :y, :z])
        @test_throws "cannot combine NamedRowArrays with different row names" hcat(an, other)
        @test_throws "cannot combine NamedRowArrays with different row names" hcat(v, an, other)
    end

    @testset "Matrix products" begin
        an = NamedRowArray(a, rownames)
        R = rand(n, n)

        # the rows of a product are the rows of the left factor
        @test R * an == R * a
        @test R * an isa Matrix
        Rn = NamedRowArray(R, rownames)
        @test Rn * an == R * a
        @test Rn * an isa NamedRowMatrix
        @test (Rn * an).rownames == an.rownames
        Rn2 = NamedRowArray(R, [:w, :x, :y, :z])
        @test (Rn2 * an).rownames == [:w, :x, :y, :z]

        x = rand(m)
        @test an * x == a * x
        @test an * x isa NamedRowVector
        @test (an * x).rownames == an.rownames
        @test (an * a').rownames == an.rownames
        @test (an * rand(m, 2)).rownames == an.rownames
        @test (an * NamedRowArray(rand(m, 2), Symbol.('a':'e'))).rownames == an.rownames
        @test (an * NamedRowArray(x, Symbol.('a':'e'))).rownames == an.rownames

        # without names on the left factor, the product is a plain array
        vn = NamedRowArray(v, rownames)
        @test R * vn == R * v
        @test R * vn isa Vector
        @test R' * an == R' * a
        @test R' * an isa Matrix
        @test v' * an == v' * a
        @test transpose(v) * an == transpose(v) * a
        @test v' * vn == v' * v
        @test transpose(v) * vn == transpose(v) * v
        @test Diagonal(v) * vn == Diagonal(v) * v
        @test Diagonal(v) * vn isa Vector
    end

    @testset "Reordering" begin
        vn = NamedRowArray([3.0, 1.0, 2.0, 4.0], rownames)
        an = NamedRowArray([3.0 2; 1 4; 2 6; 1 1], rownames)

        # the row names move with their rows
        function rows_match(B, A)
            for name in B.rownames
                B[name] == A[name] || return false
            end
            return true
        end

        @testset "reverse" begin
            @test reverse(vn) == reverse(parent(vn))
            @test reverse(vn).rownames == [:d, :c, :b, :a]
            @test reverse(vn; dims=1).rownames == [:d, :c, :b, :a]
            @test reverse(vn, 2, 3) == reverse(parent(vn), 2, 3)
            @test reverse(vn, 2, 3).rownames == [:a, :c, :b, :d]
            @test reverse(vn, 3).rownames == [:a, :b, :d, :c]
            @test rows_match(reverse(vn, 2), vn)

            @test reverse(an) == reverse(parent(an))
            @test reverse(an).rownames == [:d, :c, :b, :a]
            @test reverse(an; dims=1) == reverse(parent(an); dims=1)
            @test rows_match(reverse(an; dims=1), an)
            @test reverse(an; dims=2) == reverse(parent(an); dims=2)
            @test reverse(an; dims=2).rownames == an.rownames
        end

        @testset "sort" begin
            @test sort(vn) == sort(parent(vn))
            @test sort(vn).rownames == [:b, :c, :a, :d]
            @test sort(vn; rev=true).rownames == [:d, :a, :c, :b]
            @test rows_match(sort(vn; by=x -> -x), vn)

            @test sort(an; dims=2) == sort(parent(an); dims=2)
            @test sort(an; dims=2).rownames == an.rownames
            msg = "sorting a NamedRowMatrix along dims=1 sorts each column separately"
            @test_throws msg sort(an; dims=1)

            @test sortslices(an; dims=1) == sortslices(parent(an); dims=1)
            @test sortslices(an; dims=1).rownames == [:d, :b, :c, :a]
            @test rows_match(sortslices(an; dims=1, rev=true), an)
            @test sortslices(an; dims=2) == sortslices(parent(an); dims=2)
            @test sortslices(an; dims=2).rownames == an.rownames
            @test_throws "dims must be 1 or 2 for a NamedRowMatrix, got 3" sortslices(an; dims=3)
        end

        @testset "circshift" begin
            @test circshift(vn, 1) == circshift(parent(vn), 1)
            @test circshift(vn, 1).rownames == [:d, :a, :b, :c]
            @test circshift(vn, -5).rownames == [:b, :c, :d, :a]
            @test circshift(vn, [1]).rownames == [:d, :a, :b, :c]
            @test circshift(an, (1, 1)) == circshift(parent(an), (1, 1))
            @test circshift(an, (1, 1)).rownames == [:d, :a, :b, :c]
            @test circshift(an, (0, 1)).rownames == an.rownames
            @test rows_match(circshift(an, 3), an)
        end

        @testset "shuffle" begin
            s = shuffle(Xoshiro(1), vn)
            @test sort(collect(s)) == sort(parent(vn))
            @test rows_match(s, vn)
            @test rows_match(shuffle(vn), vn)
            msg = "shuffling a NamedRowMatrix moves elements between rows"
            @test_throws msg shuffle(an)
            @test_throws "shuffle! cannot reorder the row names" shuffle!(copy(vn))
            @test_throws "shuffle! cannot reorder the row names" shuffle!(Xoshiro(1), copy(an))
        end

        @testset "circshift!" begin
            msg = "circshift! cannot reorder the row names of a NamedRowArray in place; use circshift instead"
            @test_throws msg circshift!(copy(vn), 1)
            @test circshift!(copy(vn), 4) == vn
            @test_throws msg circshift!(similar(an), an, 1)
            @test_throws msg circshift!(similar(an), an, (1, 0))
            @test_throws msg circshift!(similar(an), an, [1, 0])
            @test circshift!(similar(an), an, (0, 1)) == circshift(parent(an), (0, 1))
            @test circshift!(similar(an), an, (4, 1)) == circshift(parent(an), (0, 1))
            @test circshift!(similar(an), an, ()) == an
            @test circshift!(similar(an), parent(an), (0, 1)).rownames == an.rownames
        end

        @testset "in place" begin
            msg = "cannot reorder the row names of a NamedRowArray in place"
            @test_throws "sort! $msg; use sort instead" sort!(copy(vn))
            @test_throws "reverse! $msg; use reverse instead" reverse!(copy(vn))
            @test_throws msg reverse!(copy(vn), 2, 3)
            @test_throws msg reverse!(copy(an))
            @test_throws msg reverse!(copy(an); dims=1)
            @test_throws "sort! $msg; use sortslices(A; dims=1) instead" sort!(copy(an); dims=1)
            @test_throws "permute! $msg; use A[p] instead" permute!(copy(vn), [2, 1, 3, 4])
            @test_throws "invpermute! $msg; use A[invperm(p)] instead" invpermute!(copy(vn), [2, 1, 3, 4])

            # reordering within rows keeps the row names valid
            anc = copy(an)
            @test reverse!(anc; dims=2) === anc
            @test anc == reverse(parent(an); dims=2)
            anc = copy(an)
            @test sort!(anc; dims=2) === anc
            @test anc == sort(parent(an); dims=2)
        end
    end

    @testset "Display" begin
        vn = NamedRowArray(v, rownames)
        an = NamedRowArray(a, rownames)
        @test summary(vn) == "4-element NamedRowVector(::Vector{Float64}) with eltype Float64"
        @test summary(an) == "4×5 NamedRowMatrix(::Matrix{Float64}) with eltype Float64"
        @test summary(view(an, [1, 3], 2:3)) ==
              "2×2 NamedRowMatrix(view(::Matrix{Float64}, [1, 3], 2:3)) with eltype Float64"
        @test summary(NamedRowArray(OffsetArray(v, 0:3), rownames)) ==
              "4-element NamedRowVector(OffsetArray(::Vector{Float64}, 0:3)) with eltype Float64 with indices 0:3"
        @test summary(transpose(an)) == "5×4 transpose(NamedRowMatrix(::Matrix{Float64})) with eltype Float64"

        # the row names are on their own line, followed by the data as Base shows them
        datalines(A) = split(sprint(show, MIME"text/plain"(), A), '\n')[2:end]
        shown = split(sprint(show, MIME"text/plain"(), an), '\n')
        @test shown[1] == "4×5 NamedRowMatrix(::Matrix{Float64}) with eltype Float64 with row names:"
        @test shown[2] == " a b c d"
        @test shown[3:end] == datalines(a)
        @test split(sprint(show, MIME"text/plain"(), vn), '\n')[3:end] == datalines(v)
        @test sprint(show, MIME"text/plain"(), NamedRowArray(Float64[], Symbol[])) ==
              "0-element NamedRowVector(::Vector{Float64}) with eltype Float64 with row names:\n "

        # names that do not fit are cut off, also within a name wider than the display
        show30(A) = sprint(show, MIME"text/plain"(), A; context=:displaysize => (24, 30))
        long = NamedRowArray(ones(3), [Symbol("x"^50), :b, :c])
        @test split(show30(long), '\n')[2] == " $("x"^28)…"
        many = NamedRowArray(ones(9), [Symbol("αβγ", i) for i in 1:9])
        @test split(show30(many), '\n')[2] == " αβγ1 αβγ2 αβγ3 αβγ4 αβγ5 …"
    end

end
