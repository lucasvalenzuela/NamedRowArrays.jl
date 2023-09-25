using NamedRowArrays
using Test

using FillArrays

fill4(v::AbstractVector, idx) = fill(4, count(!=(false), idx))
fill4(v::AbstractVector, idx::Colon) = fill(4, length(v))

fill4(a::AbstractMatrix, idx, idx2) = fill(4, count(!=(false), idx), count(!=(false), idx2))
fill4(a::AbstractMatrix, idx::Colon, idx2) = fill(4, size(a, 1), count(!=(false), idx2))
fill4(a::AbstractMatrix, idx, idx2::Colon) = fill(4, count(!=(false), idx), size(a, 2))
fill4(a::AbstractMatrix, idx::Colon, idx2::Colon) = fill(4, size(a, 1), size(a, 2))

@testset "NamedRowArrays.jl" begin
    n = 4
    m = 4
    v = rand(n)
    a = rand(n, m)
    rownames = Symbol.('a':'z')[1:n]

    arrs = [
          1:n,
          rand(m, n)',
          transpose(rand(m, n)),
          Fill(1.0, 4, m),
          @view rand(2n, m)[1:n, :]
         ]

    idxs = [
            (:),
            1:n,
            2:3,
            3:3,
            n:-1:1,
            [1,3],
            [3,1],
            [3],
            isodd.(1:n),
            collect(isodd.(1:n)),
           ]

    idxs2 = [
            (:),
            1:m,
            2:3,
            3:3,
            m:-1:1,
            [1,3],
            [3,1],
            [3],
            isodd.(1:m),
            collect(isodd.(1:m)),
           ]

    @testset "Creating" begin
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

        @test_throws ArgumentError NamedRowArray(v, rownames[1:3])
        @test_throws ArgumentError NamedRowArray(a, rownames[1:3])
        @test_throws ArgumentError NamedRowArray(v, String.(rownames))
        @test_throws ArgumentError NamedRowArray(a, String.(rownames))
        @test_throws ArgumentError NamedRowArray(rand(4, 4, 4), rownames)

        @test_throws MethodError NamedRowArray(1, rownames)
        @test_throws MethodError NamedRowVector(a, rownames)
        @test_throws MethodError NamedRowMatrix(v, rownames)

        @test names(vn) == vn.rownames
        @test names(an) == an.rownames

        for arr in arrs
            @test_nowarn NamedRowArray(arr, rownames)
        end
    end

    @testset "Array Functions" begin
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

        @test copy(an) !== an
        @test copy(an).data !== an.data
        @test copy(an).rownames == an.rownames
        @test deepcopy(an) !== an
    end

    @testset "Indexing" begin
        vn = NamedRowArray(v, rownames)
        an = NamedRowArray(a, rownames)

        @testset "Vector" begin
            # single elements
            @test vn[1] == v[1]
            @test vn[2] == v[2]
            @test_throws BoundsError vn[0]
            @test_throws BoundsError vn[n+1]

            @test vn[rownames[1]] == v[1]
            @test vn[rownames[2]] == v[2]
            @test_throws ArgumentError vn[:foo]
            @test_throws ArgumentError vn["foo"]

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
            @test_throws ArgumentError vn[[:foo, :a]]
            @test_throws ArgumentError vn[["foo", "a"]]
        end

        @testset "Matrix" begin
            # single elements
            @test an[1] == a[1]
            @test an[2] == a[2]
            @test an[n+1] == a[n+1]
            @test_throws BoundsError an[0]
            @test_throws BoundsError an[length(an)+1]

            @test an[1,1] == a[1,1]
            @test an[2,3] == a[2,3]
            @test an[n,3] == a[n,3]
            @test_throws BoundsError an[n+1,1]
            @test_throws BoundsError an[0,1]

            @test an[CartesianIndex(1,2)] == a[CartesianIndex(1,2)]
            @test an[[CartesianIndex(1,2),CartesianIndex(2,2)]] == a[[CartesianIndex(1,2),CartesianIndex(2,2)]]

            @test an[rownames[1],1] == a[1,1]
            @test an[rownames[2],3] == a[2,3]
            @test_throws ArgumentError an[:foo,1]
            @test_throws ArgumentError an["foo",1]

            # multiple elements
            for idx in idxs
                for idx2 in idxs2
                    @test an[idx,idx2] == a[idx,idx2]
                    @test an[idx,idx2] isa NamedRowMatrix
                    @test all(an[idx,idx2].rownames .== rownames[idx])
                end

                @test an[idx,2] == a[idx,2]
                @test an[idx,2] isa NamedRowVector
                @test all(an[idx,2].rownames .== rownames[idx])
            end

            for idx2 in idxs2
                @test an[rownames[1],idx2] == a[1,idx2]
                @test an[rownames[2],idx2] == a[2,idx2]
                @test_throws ArgumentError an[:foo,idx2]
                @test_throws ArgumentError an["foo",idx2]

                @test an[rownames[2:3],idx2] == a[2:3,idx2]
                @test an[rownames[2:3],idx2] isa NamedRowMatrix
                @test an[rownames[2:3],idx2].rownames == an.rownames[2:3]
                @test_throws ArgumentError an[[:foo, :a],idx2]
                @test_throws ArgumentError an[["foo", "a"],idx2]
            end

            @test an[rownames[1]] == a[1,:]
            @test an[rownames[2]] == a[2,:]
            @test an[rownames[1]] isa Vector
            @test_throws ArgumentError an[:foo]
            @test_throws ArgumentError an["foo"]

            @test an[rownames[2:3],1] == a[2:3,1]
            @test an[rownames[2:3],3] == a[2:3,3]
            @test an[rownames[2:3],1] isa NamedRowVector
            @test an[rownames[2:3],1].rownames == an.rownames[2:3]
            @test_throws ArgumentError an[[:foo, :a],1]
            @test_throws ArgumentError an[["foo", "a"],1]
        end
    end

    @testset "Setting" begin
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
            @test_throws BoundsError (vnc[n+1] = 4)

            vnc = copy(vn)
            vnc[rownames[1]] = 4
            @test vnc[rownames[1]] == 4
            vnc[rownames[2]] = 4
            @test vnc[rownames[2]] == 4
            @test_throws ArgumentError vnc[:foo]
            @test_throws ArgumentError vnc["foo"]

            # multiple elements
            for idx in idxs
                vnc = copy(vn)
                vnc[idx] = fill4(vnc, idx)
                @test all(vnc[idx] .== 4)

                vnc = copy(vn)
                vnc[rownames[2:3]] = fill4(vnc, 2:3)
                @test all(vnc[rownames[2:3]] .== 4)
                @test_throws ArgumentError (vnc[[:foo, :a]] = fill4(vnc, 2:3))
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
            anc[n+1] = 4
            @test anc[n+1] == 4
            @test_throws BoundsError (anc[0] = 4)
            @test_throws BoundsError (anc[length(an)+1] = 4)

            anc = copy(an)
            anc[1,1] = 4
            @test anc[1,1] == 4
            anc[2,3] = 4
            @test anc[2,3] == 4
            anc[n,3] = 4
            @test anc[n,3] == 4
            anc[CartesianIndex(3,2)] = 4
            @test anc[CartesianIndex(3,2)] == 4
            @test_throws BoundsError (anc[n+1,1] = 4)
            @test_throws BoundsError (anc[0,1] = 4)

            anc = copy(an)
            anc[rownames[1],1] = 4
            @test anc[rownames[1],1] == 4
            anc[rownames[2],3] = 4
            @test anc[rownames[2],3] == 4
            @test_throws ArgumentError (anc[:foo,1] = 4)
            @test_throws ArgumentError (anc["foo",1] = 4)

            # multiple elements
            for idx in idxs
                for idx2 in idxs2
                    anc = copy(an)
                    anc[idx,idx2] = fill4(anc, idx, idx2)
                    @test all(anc[idx,idx2] .== 4)
                end

                anc = copy(an)
                anc[idx,2] = fill4(anc, idx, 1)
                @test all(anc[idx,2] .== 4)
            end

            for idx2 in idxs2
                anc = copy(an)
                anc[rownames[1],idx2] = fill4(anc, 1, idx2)
                @test all(anc[rownames[1],idx2] .== 4)
                anc[rownames[2],idx2] = fill4(anc, 1, idx2)
                @test all(anc[rownames[2],idx2] .== 4)
                @test_throws ArgumentError (anc[:foo,idx2] = fill4(anc, 1, idx2))
                @test_throws ArgumentError (anc["foo",idx2] = fill4(anc, 1, idx2))

                anc = copy(an)
                anc[rownames[2:3],idx2] = fill4(anc, 2:3, idx2)
                @test all(anc[rownames[2:3],idx2] .== 4)
                @test_throws ArgumentError (anc[[:foo, :a],idx2] = fill4(anc, 2:3, idx2))
                @test_throws ArgumentError (anc[["foo", "a"],idx2] = fill4(anc, 2:3, idx2))
            end

            anc = copy(an)
            anc[[CartesianIndex(3,2), CartesianIndex(3,3)]] = [4, 4]
            @test all(anc[[CartesianIndex(3,2), CartesianIndex(3,3)]] .== 4)

            anc = copy(an)
            anc[rownames[1]] = fill4(anc, 1, :)
            @test all(anc[rownames[1]] .== 4)
            anc[rownames[2]] = fill4(anc, 1, :)
            @test all(anc[rownames[2]] .== 4)
            @test_throws ArgumentError (anc[:foo] = 4)
            @test_throws ArgumentError (anc["foo"] = 4)

            anc = copy(an)
            anc[rownames[2:3]] = fill4(anc, 2:3, :)
            @test all(anc[rownames[2:3]] .== 4)
            @test_throws ArgumentError (anc[:foo, :a] = fill4(anc, 2:3, :))
            @test_throws ArgumentError (anc["foo", "a"] = fill4(anc, 2:3, :))
        end
    end
end
