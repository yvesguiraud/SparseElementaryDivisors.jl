using Test, Random, SparseArrays
using SparseElementaryDivisors
import Hecke

# Reference: Hecke on the dense matrix, zeros removed, increasing order.
function reference(A::AbstractMatrix{<:Integer})
    M = Hecke.matrix(Hecke.ZZ, Matrix{Int}(A))
    return sort!(filter(!iszero, BigInt.(Hecke.elementary_divisors(M))))
end

# Random matrix (proportion `p_unit` of ±1 entries, the others in ±2:max_abs).
function random_matrix(m, n, density; p_unit = 0.6, max_abs = 5, seed = 1)
    rng = MersenneTwister(seed)
    A = zeros(Int, m, n)
    for r ∈ 1:m, c ∈ 1:n
        if rand(rng) < density
            v = rand(rng) < p_unit ? 1 : rand(rng, 2:max_abs)
            rand(rng) < 0.5 && (v = -v)
            A[r, c] = v
        end
    end
    return A
end

# Same generation, row by row (`rows` interface).
function random_rows_data(m, n, density; p_unit = 0.6, max_abs = 5, seed = 1)
    rng = MersenneTwister(seed)
    rows = Vector{Vector{Tuple{Int,Int}}}(undef, m)
    for r ∈ 1:m
        row = Tuple{Int,Int}[]
        for c ∈ 1:n
            if rand(rng) < density
                v = rand(rng) < p_unit ? 1 : rand(rng, 2:max_abs)
                rand(rng) < 0.5 && (v = -v)
                push!(row, (c, v))
            end
        end
        rows[r] = row
    end
    return rows
end

const NW = min(4, Threads.nthreads())
const sed = sparse_elementary_divisors

# True if `f()` throws an ArgumentError whose message contains `needle` (`@test_throws`
# with a string checks the message only, not the type).
function throws_argument_error(f, needle)
    err = try
        f()
        nothing
    catch e
        e
    end
    return err isa ArgumentError && occursin(needle, sprint(showerror, err))
end

@testset "SparseElementaryDivisors" begin

    @testset "simple cases" begin
        @test sed(sparse(zeros(Int, 3, 4))) == BigInt[]
        @test sed(sparse([1 0 0; 0 1 0; 0 0 1])) == BigInt[1, 1, 1]
        @test sed(sparse([2 0; 0 3])) == BigInt[1, 6]
        @test sed(sparse([2 4; 6 8])) == reference([2 4; 6 8])
    end

    @testset "random matrices against Hecke" begin
        for (m, n, dens, seed) ∈ ((30, 40, 0.15, 1), (60, 45, 0.10, 2), (50, 50, 0.20, 3),
                                   (40, 80, 0.08, 4), (80, 40, 0.08, 5))
            A = random_matrix(m, n, dens; seed)
            ref = reference(A)
            @test sed(sparse(A); value_type = Int64) == ref
            @test sed(sparse(A); value_type = Int64, nworkers = 1) == ref
            @test sed(sparse(A); value_type = Int64, nworkers = NW) == ref
        end
    end

    @testset "core and pivoting options" begin
        A = random_matrix(60, 80, 0.10; seed = 7)
        ref = reference(A)
        for fmt ∈ (:auto, :sparse, :dense)
            @test sed(sparse(A); value_type = Int64, core_format = fmt) == ref
        end
        for b ∈ (1, 2, typemax(Int))
            @test sed(sparse(A); value_type = Int64, max_buckets_scanned = b) == ref
        end
    end

    @testset "input types" begin
        A = random_matrix(40, 50, 0.12; seed = 11)
        ref = reference(A)
        rows = random_rows_data(40, 50, 0.12; seed = 11)
        @test sed(rows, 40, 50; value_type = Int64) == ref
        @test sed(40, 50, rows; value_type = Int64) == ref
        @test sed(A; value_type = Int64) == ref                                   # dense
        @test sed(A'; value_type = Int64) == reference(Matrix(A'))                # Adjoint
        @test sed(view(sparse(A), 1:30, 1:40); value_type = Int64) == reference(A[1:30, 1:40])
        S = Hecke.sparse_matrix(Hecke.ZZ)                                         # SMat de Hecke
        for r ∈ 1:40
            push!(S, Hecke.sparse_row(Hecke.ZZ, [(c, Hecke.ZZ(A[r, c])) for c ∈ 1:50 if A[r, c] != 0]))
        end
        @test sed(S; value_type = Int64) == ref
        @test sed(Hecke.matrix(Hecke.ZZ, A); value_type = Int64) == ref          # ZZMatrix
    end

    @testset "invalid input" begin
        R = Vector{Tuple{Int,Int}}
        bad = throws_argument_error
        # a valid 2 × 3 input, for comparison
        @test sed([[(1, 1), (3, 2)], [(2, -1)]], 2, 3) == reference([1 0 2; 0 -1 0])
        # columns out of range
        @test bad(() -> sed([[(1, 1), (0, 2)], [(2, -1)]], 2, 3), "row 1 (entry 2): column 0 is outside 1:3")
        @test bad(() -> sed([[(1, 1)], [(4, 1)]], 2, 3), "row 2 (entry 1): column 4 is outside 1:3")
        @test bad(() -> sed([[(-2, 1)], [(2, 1)]], 2, 3), "row 1 (entry 1): column -2 is outside 1:3")
        # unsorted and repeated columns
        @test bad(() -> sed([[(3, 1), (1, 2)], [(2, 1)]], 2, 3), "row 1 (entry 2): column 1 after column 3")
        @test bad(() -> sed([[(1, 1)], [(2, 1), (2, 5)]], 2, 3), "row 2 (entry 2): column 2 appears twice")
        @test bad(() -> sed([[(1, 1)], [(2, 1), (2, 5)]], 2, 3), "sparse(I, J, V) sorts and sums duplicates")
        # size
        @test bad(() -> sed([[(1, 1)], [(2, 1)], [(3, 1)]], 2, 3), "rows has 3 rows but m = 2")
        @test bad(() -> sed(R[], -1, 3), "m and n must be nonnegative (got m = -1, n = 3)")
        @test bad(() -> sed([R()], 1, -1), "m and n must be nonnegative (got m = 1, n = -1)")
        # the form (m, n, rows) is checked too
        @test bad(() -> sed(2, 3, [[(1, 1)], [(4, 1)]]), "row 2 (entry 1): column 4 is outside 1:3")
        # a non-canonical SparseMatrixCSC (row 1 twice in column 1) gives a repeated column
        A = SparseMatrixCSC(2, 2, [1, 3, 4], [1, 1, 2], [1, 1, 1])
        @test bad(() -> sed(A), "row 1 (entry 2): column 1 appears twice")
        # valid limit cases: empty rows, no row, no column, zero values (ignored)
        @test sed([R(), R()], 2, 3) == BigInt[]
        @test sed(R[], 0, 3) == BigInt[]
        @test sed([R(), R()], 2, 0) == BigInt[]
        @test sed([[(1, 0), (2, 3)], [(1, 2), (3, 0)]], 2, 3) == reference([0 3 0; 2 0 0])
    end

    @testset "Matrix Market" begin
        A = random_matrix(30, 35, 0.15; seed = 13)
        ref = reference(A)
        coord = IOBuffer()
        println(coord, "%%MatrixMarket matrix coordinate integer general")
        println(coord, "% un commentaire")
        I, J, V = findnz(sparse(A))
        println(coord, "30 35 ", length(V))
        for k ∈ eachindex(V)
            println(coord, I[k], " ", J[k], " ", V[k])
        end
        f = tempname() * ".mtx"
        write(f, take!(coord))
        @test read_matrix_market(f) == sparse(A)
        @test sed(f; value_type = Int64) == ref
        rm(f)
        # symmetric: only the lower triangle is stored
        io = IOBuffer("%%MatrixMarket matrix coordinate integer symmetric\n3 3 4\n1 1 2\n2 1 1\n3 2 3\n3 3 5\n")
        @test read_matrix_market(io) == sparse([2 1 0; 1 0 3; 0 3 5])
        # array format, column by column
        io = IOBuffer("%%MatrixMarket matrix array integer general\n2 3\n1\n0\n0\n4\n5\n0\n")
        @test read_matrix_market(io) == sparse([1 0 5; 0 4 0])
        @test_throws ArgumentError read_matrix_market(IOBuffer("pas un fichier Matrix Market\n"))
    end

    @testset "explicit overflow" begin
        rows = random_rows_data(120, 120, 0.15; seed = 1)
        @test_throws ErrorException sed(rows, 120, 120; value_type = Int8, nworkers = NW)
        # after an overflow, a rerun with Int64 succeeds
        @test sed(rows, 120, 120; value_type = Int64, nworkers = NW) isa Vector{BigInt}
    end

    @testset "entries out of range" begin
        @test_throws ArgumentError sed(sparse([big(2)^80 0; 0 1]))
    end

    @testset "Magma (script generation and output parsing)" begin
        rows = [[(1, 2)], [(2, 3), (3, -1)]]
        @test SparseElementaryDivisors._magma_matrix_definition(rows, 2, 3) ==
              "M := SparseMatrix(2, 3, [ 1, 1, 2, 2, 2, 3, 3, -1 ]);"
        raw = "Time: 0.630\n##DIVS_BEGIN##\n[ 1, 1, 6,\n 12345678901234567890 ]\n##DIVS_END##\n"
        divs, t = SparseElementaryDivisors._parse_magma_output(raw)
        @test divs == BigInt[1, 1, 6, big"12345678901234567890"]
        @test t == 0.630
        @test_throws ErrorException SparseElementaryDivisors._parse_magma_output(
            "##MAGMA_ERROR_BEGIN##Runtime error##MAGMA_ERROR_END##")
    end
end

include("test_fingerprints.jl")