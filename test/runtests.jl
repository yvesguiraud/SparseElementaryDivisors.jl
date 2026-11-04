using Test, Random, SparseArrays
using SparseElementaryDivisors
import Hecke
using BitIntegers: Int256

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

# A sparse matrix type other than SparseMatrixCSC (as other packages define), for the
# conversion path of sparse_elementary_divisors.
struct WrappedSparse <: AbstractSparseMatrix{Int, Int}
    A::SparseMatrixCSC{Int, Int}
end
Base.size(W::WrappedSparse) = size(W.A)
Base.getindex(W::WrappedSparse, i::Int, j::Int) = W.A[i, j]

const NW = min(4, Threads.nthreads())
const sed = sparse_elementary_divisors

# The exception thrown by `f()`, or `nothing`.
function caught(f)
    try
        f()
    catch e
        return e
    end
    return nothing
end

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
            @test sed(sparse(A); value_type = BigInt) == ref
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
        @test sed(WrappedSparse(sparse(A)); value_type = Int64) == ref            # other sparse type
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
        println(coord, "% a comment")
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
        @test_throws ArgumentError read_matrix_market(IOBuffer("not a Matrix Market file\n"))
        # real field: integer values only
        io = IOBuffer("%%MatrixMarket matrix coordinate real general\n2 2 2\n1 1 2.0\n2 2 -3\n")
        @test read_matrix_market(io) == sparse([2 0; 0 -3])
        io = IOBuffer("%%MatrixMarket matrix coordinate real general\n2 2 1\n1 1 2.5\n")
        @test throws_argument_error(() -> read_matrix_market(io), "non-integer value \"2.5\"")
        io = IOBuffer("%%MatrixMarket matrix vector integer general\n2 2\n")
        @test throws_argument_error(() -> read_matrix_market(io), "unsupported format \"vector\"")
    end

    @testset "Matrix Market writing" begin
        header = "%%MatrixMarket matrix coordinate integer general\n"
        # the 24 files of test/data, read then written with their comment: the same bytes
        data = joinpath(@__DIR__, "data")
        files = sort!(filter(endswith(".mtx"), readdir(data)))
        @test length(files) == 24
        for f ∈ files
            bytes = read(joinpath(data, f))
            comment = [l[3:end] for l ∈ split(String(copy(bytes)), '\n') if startswith(l, "% ")]
            io = IOBuffer()
            write_matrix_market(io, read_matrix_market(joinpath(data, f)); comment)
            @test take!(io) == bytes
        end
        # dense input, by path, comment of two lines; read back
        A = random_matrix(30, 35, 0.15; seed = 13)
        path = tempname() * ".mtx"
        @test write_matrix_market(path, A; comment = "first line\nsecond line\n") === nothing
        @test read_matrix_market(path) == sparse(A)
        text = read(path, String)
        @test startswith(text, header * "% first line\n% second line\n30 35 $(count(!iszero, A))\n")
        @test endswith(text, "\n") && !occursin('\r', text)
        rm(path)
        # stored zero left out, entries by row then column, vector of comments
        S = SparseMatrixCSC(2, 2, [1, 3, 4], [1, 2, 1], [-3, 0, 5])
        @test sprint(write_matrix_market, S) == header * "2 2 2\n1 1 -3\n1 2 5\n"
        @test sprint(io -> write_matrix_market(io, S; comment = ["a", "b\r\nc"])) ==
              header * "% a\n% b\n% c\n2 2 2\n1 1 -3\n1 2 5\n"
        # no entry; values beyond a machine word; Bool
        E = spzeros(Int, 3, 0)
        @test sprint(write_matrix_market, E) == header * "3 0 0\n"
        @test read_matrix_market(IOBuffer(sprint(write_matrix_market, E))) == E
        @test sprint(write_matrix_market, sparse([big(2)^70 0; 0 -1])) ==
              header * "2 2 2\n1 1 1180591620717411303424\n2 2 -1\n"
        @test sprint(write_matrix_market, sparse([true false; false true])) == header * "2 2 2\n1 1 1\n2 2 1\n"
        @test_throws ArgumentError write_matrix_market(IOBuffer(), S; comment = 3)
    end

    @testset "degree buckets" begin
        SED = SparseElementaryDivisors
        # largest degree with its own bucket: clamp(max(m, n), 1, 10_000)
        for (m, n, expected) ∈ ((1, 56, 56), (0, 0, 1), (1, 10_001, 10_000), (300, 200, 300))
            rows = [Tuple{Int,Int}[] for _ ∈ 1:m]
            M = SED.DualSparseMatrix(rows, m, n; value_type = Int64)
            @test M.row_buckets.max_deg == M.col_buckets.max_deg == M.unit_index.buckets.max_deg == expected
            @test length(M.unit_index.buckets.buckets) == expected + 2
            @test SED.check_structure_consistency(M)
        end
        # the cursor low and the rest of the structure stay consistent along a reduction
        checks = Bool[]
        for (m, n, dens, seed) ∈ ((60, 80, 0.10, 7), (80, 40, 0.08, 5), (50, 50, 0.20, 3))
            A = random_matrix(m, n, dens; seed)
            M = SED.DualSparseMatrix(SED._rows_data(sparse(A)), m, n; value_type = Int64)
            probe = () -> push!(checks, SED.check_structure_consistency(M))
            SED.reduce_to_core!(M; probe)
            @test SED.check_structure_consistency(M)
        end
        @test !isempty(checks) && all(checks)
    end

    @testset "internals" begin
        SED = SparseElementaryDivisors
        A = random_matrix(60, 80, 0.10; seed = 7)
        ref = reference(A)
        rows = SED._rows_data(sparse(A))
        # overflow bucket: with max_tracked_degree = 3, the degrees above 3 share it
        M = SED.DualSparseMatrix(rows, 60, 80; value_type = Int64, max_tracked_degree = 3)
        @test SED._reduce_and_divisors!(M) == ref
        @test SED.check_structure_consistency(M)
        # display of the trajectory fingerprint
        tr = SED.ReductionTrace()
        sed(sparse(A); trace = tr)
        @test startswith(sprint(show, tr), "ReductionTrace(h=0x") && occursin("pivots=$(tr.n_pivots)", sprint(show, tr))
        # sparse extraction of an empty core
        M = SED.DualSparseMatrix(SED._rows_data(sparse([1 0; 0 1])), 2, 2; value_type = Int64)
        @test SED.reduce_to_core!(M) == 2
        @test Hecke.nrows(SED.extract_core_smat(M)) == 0
        # unwrapping of a TaskFailedException
        t = Threads.@spawn error("inner failure")
        e = caught(() -> wait(t))
        @test e isa TaskFailedException
        @test SED._unwrap_task_failure(e) isa ErrorException
    end

    @testset "diagnostics" begin
        SED = SparseElementaryDivisors
        A = random_matrix(40, 50, 0.12; seed = 11)
        ref = reference(A)
        rows = SED._rows_data(sparse(A))
        quiet(f) = redirect_stdout(f, devnull)
        # at each step, the diagnostic search (full loop over the buckets) finds the pivot of
        # the search of the reduction (cursor of P5), overflow bucket included
        for cfg ∈ (ReductionConfig(nworkers = 1), ReductionConfig(nworkers = 1, max_buckets_scanned = 3)),
            mtd ∈ (SED._max_tracked_degree(40, 50), 3)
            M = SED.DualSparseMatrix(rows, 40, 50; value_type = Int64, max_tracked_degree = mtd)
            wsp = SED.ReductionWorkspace{Int64}(cfg)
            same = true
            steps = 0
            while true
                SED.eliminate_structured!(M)
                d = SED.diagnostic_pivot_search(M; config = cfg)
                p_r, p_c, score = SED._search_pivot!(M, wsp, cfg, 1, nothing)
                same &= (d.p_r, d.p_c, d.best_score) == (p_r, p_c, score)
                d.found || break
                SED.markowitz_pivot!(M; config = cfg, workspace = wsp) || break
                steps += 1
            end
            @test same && steps > 0
            @test SED.check_structure_consistency(M)
        end
        # cost of the searches, its summary, sweep of two settings (printed output discarded)
        M = SED.DualSparseMatrix(rows, 40, 50; value_type = Int64)
        diag = SED.diagnose_pivot_search_cost(deepcopy(M))
        @test !isempty(diag)
        totals = quiet(() -> SED.summarize_pivot_search_diag(diag; top_k = 3))
        @test totals.n_total == length(diag)
        results = quiet(() -> SED.sweep_search_params(M, [(4, 1), (0, typemax(Int))]))
        @test length(results) == 2 && all(r -> r.divisors == ref, results)
        I2 = SED.DualSparseMatrix(SED._rows_data(sparse([1 0; 0 1])), 2, 2; value_type = Int64)
        @test only(quiet(() -> SED.sweep_search_params(I2, [(4, 1)]))).divisors == [1, 1]  # empty core
        # the checks detect a corrupted structure
        @test SED.check_structure_consistency(M)
        unit_col = findfirst(c -> M.unit_index.unit_count[c] > 0 && length(M.col_rows[c]) >= 2, 1:50)
        row = findfirst(r -> !isempty(M.row_cols[r]), 1:40)
        corruptions = (C -> push!(C.row_vals[row], 1),                                # row lengths
                       C -> (deleteat!(C.col_rows[unit_col], 1); deleteat!(C.col_vals[unit_col], 1)),
                       C -> (C.unit_index.unit_count[unit_col] += 1),                  # unit count
                       C -> (C.unit_index.unit_count[unit_col] = 0),                   # tracked, no unit
                       C -> SED._bucket_remove!(C.unit_index.buckets, unit_col),       # untracked
                       C -> (C.row_buckets.pos_in_bucket[row] = 0),                    # position
                       C -> (C.unit_index.buckets.low = C.unit_index.buckets.max_deg + 1))
        for corrupt! ∈ corruptions
            C = deepcopy(M)
            corrupt!(C)
            @test !quiet(() -> SED.check_structure_consistency(C))
        end
    end

    @testset "explicit overflow" begin
        rows = random_rows_data(120, 120, 0.15; seed = 1)
        @test_throws CoefficientOverflowError sed(rows, 120, 120; value_type = Int8, nworkers = NW)
        err = caught(() -> sed(rows, 120, 120; value_type = Int8, nworkers = NW))
        @test err.value_type === Int8 && err.entry === nothing && err.unit_pivots > 0
        @test occursin("Arithmetic overflow detected with value_type=Int8, after $(err.unit_pivots) unit pivots",
                       sprint(showerror, err))
        # after an overflow, a rerun with Int64 succeeds; the default widens by itself
        d64 = sed(rows, 120, 120; value_type = Int64, nworkers = NW)
        @test d64 isa Vector{BigInt}
        @test sed(rows, 120, 120; nworkers = NW) == d64
        @test sed(rows, 120, 120; value_type = (Int8, Int16, Int32, Int64), nworkers = NW) == d64
    end

    @testset "value types and widening" begin
        bad = throws_argument_error
        A = sparse([2 4; 6 8])
        for vt ∈ (UInt64, Float64, :fast, (), (Int64, Int32), (Int64, Int64), (BigInt, Int64), [Int256, Int128])
            @test bad(() -> sed(A; value_type = vt), "value_type")
        end
        @test sed(A; value_type = [Int32, Int64]) == reference([2 4; 6 8])
        # [1 a; a 0], a = 2^40: Int32 cannot hold a, Int64 overflows on a² = 2^80 (core −2^80)
        a = 2^40
        B = sparse([1 a; a 0])
        expected = BigInt[1, big(2)^80]
        @test sed(B) == expected                                   # :auto: Int64, then Int128
        for fmt ∈ (:auto, :dense, :sparse)
            @test sed(B; value_type = Int128, core_format = fmt) == expected
        end
        @test sed(B; value_type = (Int64, BigInt)) == expected
        @test sed(B; value_type = (Int64, Int256)) == expected     # BitIntegers.jl
        err = caught(() -> sed(B; value_type = Int32))
        @test err isa CoefficientOverflowError && err.value_type === Int32 &&
              err.unit_pivots == 0 && err.entry == (1, 2, big(2)^40)
        @test occursin("input entry (1, 2) = 1099511627776 does not fit", sprint(showerror, err))
        err = caught(() -> sed(B; value_type = Int64))
        @test err isa CoefficientOverflowError && err.value_type === Int64 && err.entry === nothing
        # Int128 with a dense core of small values
        C = random_matrix(60, 80, 0.10; seed = 7)
        @test sed(sparse(C); value_type = Int128, core_format = :dense) == reference(C)
        # verbose: the failed attempt and the type retained
        out = mktemp() do path, io
            redirect_stdout(io) do
                sed(B; verbose = true)
            end
            flush(io)
            read(path, String)
        end
        @test occursin("value_type = Int64: overflow after 0 unit pivots", out)
        @test occursin("retrying with Int128", out) && occursin("value_type = Int128 retained", out)
    end

    @testset "entries beyond a machine word" begin
        b = big(2)^70
        expected = BigInt[1, 3 * b]                                # diag(2^70, 3)
        @test sed(sparse([b 0; 0 3])) == expected
        @test sed([[(1, b)], [(2, big(3))]], 2, 2) == expected
        S = Hecke.sparse_matrix(Hecke.ZZ)
        push!(S, Hecke.sparse_row(Hecke.ZZ, [(1, Hecke.ZZ(b))]))
        push!(S, Hecke.sparse_row(Hecke.ZZ, [(2, Hecke.ZZ(3))]))
        @test sed(S) == expected
        @test sed(Hecke.matrix(Hecke.ZZ, [b 0; 0 3])) == expected
        err = caught(() -> sed(sparse([b 0; 0 3]); value_type = (Int32, Int64)))
        @test err isa CoefficientOverflowError && err.value_type === Int64 && err.entry == (1, 1, b)
        @test sed(sparse([-b 0; 0 3]); value_type = Int128) == expected
    end

    @testset "ReductionConfig" begin
        bad = throws_argument_error
        c = ReductionConfig()
        @test (c.acceptable_score_threshold, c.max_buckets_scanned, c.nworkers) == (4, 1, nothing)
        @test (c.apply_min_work, c.apply_chunk_max, c.search_mode, c.search_parallel_min_work) ==
              (25_000, 32, :auto, 5000)
        @test (c.search_ns_per_entry, c.parallel_fixed_cost_ns, c.apply_parallel_penalty,
               c.apply_ema_alpha) == (8.0, 50_000.0, 0.1, 0.05)
        @test (c.search_chunk_min, c.search_chunk_max) == (4, 256)
        @test (c.core_format, c.dense_min_density, c.dense_max_bytes, c.release_empty_vectors) ==
              (:auto, 0.05, 16 * 2^30, false)
        # display: only the fields that differ from the default, as a valid constructor call
        @test sprint(show, c) == "ReductionConfig()"
        c2 = ReductionConfig(nworkers = 4, search_mode = :always)
        @test sprint(show, c2) == "ReductionConfig(nworkers = 4, search_mode = :always)"
        @test eval(Meta.parse(sprint(show, c2))) === c2
        # checks at construction
        for (kw, needle) ∈ ((acceptable_score_threshold = -1,) => "acceptable_score_threshold must be ≥ 0",
                            (max_buckets_scanned = 0,) => "max_buckets_scanned must be ≥ 1",
                            (nworkers = 0,) => "nworkers must be nothing or ≥ 1 (got 0)",
                            (apply_min_work = -1,) => "apply_min_work must be ≥ 0",
                            (apply_chunk_max = 0,) => "apply_chunk_max must be ≥ 1",
                            (search_mode = :alway,) => "search_mode must be :auto, :always, :never or :floor (got :alway)",
                            (search_parallel_min_work = -1,) => "search_parallel_min_work must be ≥ 0",
                            (search_ns_per_entry = -1.0,) => "search_ns_per_entry must be ≥ 0",
                            (parallel_fixed_cost_ns = -1.0,) => "parallel_fixed_cost_ns must be ≥ 0",
                            (apply_parallel_penalty = -0.1,) => "apply_parallel_penalty must be ≥ 0",
                            (apply_ema_alpha = 0.0,) => "apply_ema_alpha must be in (0, 1]",
                            (apply_ema_alpha = 1.5,) => "apply_ema_alpha must be in (0, 1]",
                            (search_chunk_min = 0,) => "search_chunk_min must be ≥ 1",
                            (search_chunk_min = 8, search_chunk_max = 4) => "search_chunk_max must be ≥ search_chunk_min",
                            (core_format = :csc,) => "core_format must be :auto, :sparse or :dense",
                            (dense_min_density = 1.5,) => "dense_min_density must be in [0, 1]",
                            (dense_max_bytes = -1,) => "dense_max_bytes must be ≥ 0")
            @test bad(() -> ReductionConfig(; kw...), needle)
        end
        # keywords of sparse_elementary_divisors: a field given directly wins over config
        A = sparse(random_matrix(60, 80, 0.10; seed = 7))
        ref = reference(Matrix(A))
        st = ReductionStats()
        @test sed(A; config = ReductionConfig(nworkers = 2), stats = st) == ref
        @test st.nworkers == 2
        @test sed(A; config = ReductionConfig(nworkers = 2), nworkers = 1, stats = st) == ref
        @test st.nworkers == 1
        @test sed(A; stats = st) == ref
        @test st.nworkers == Threads.nthreads()
        @test sed(A; value_type = Int64, search_mode = :never, apply_min_work = 0, core_format = :sparse,
                  release_empty_vectors = true, acceptable_score_threshold = 0) == ref
        @test bad(() -> sed(A; nworker = 2), "unknown keyword argument nworker")
        @test bad(() -> sed(A; core_format = :csc), "core_format must be")
    end

    @testset "ReductionStats" begin
        st = ReductionStats()
        @test st.value_type === nothing && isempty(st.failed_attempts)
        @test (st.unit_pivots, st.core_size, st.core_format, st.nworkers) == (0, (rows = 0, cols = 0, nnz = 0), :none, 0)
        # empty core: two unit pivots of the structured elimination
        @test sed(sparse([1 0; 0 1]); stats = st) == BigInt[1, 1]
        @test (st.value_type, st.unit_pivots, st.markowitz_pivots) == (Int32, 2, 0)
        @test (st.core_size, st.core_format) == ((rows = 0, cols = 0, nnz = 0), :none)
        @test sprint(show, st) == "ReductionStats(value_type = Int32, unit_pivots = 2, core = 0×0)"
        plain = sprint(show, MIME"text/plain"(), st)
        @test occursin("unit pivots    2 (Markowitz 0)", plain) && occursin("core           0×0, 0 entries, none", plain)
        # [1 a; a 0], a = 2^40: one Markowitz pivot, then a 1 × 1 core −2^80
        a = 2^40
        B = sparse([1 a; a 0])
        expected = BigInt[1, big(2)^80]
        for fmt ∈ (:sparse, :dense)
            @test sed(B; value_type = Int128, core_format = fmt, stats = st) == expected
            @test (st.value_type, st.unit_pivots, st.markowitz_pivots) == (Int128, 1, 1)
            @test (st.core_size, st.core_format) == ((rows = 1, cols = 1, nnz = 1), fmt)
            @test st.markowitz_pivots == st.parallel_applications + st.sequential_applications
            @test all(>=(0), (st.sge_seconds, st.search_seconds, st.apply_seconds, st.core_seconds))
            @test isempty(st.failed_attempts)
        end
        # widening: Int32 skipped (the input does not fit), Int64 overflows, Int128 retained
        @test sed(B; stats = st) == expected
        @test st.value_type === Int128 && st.core_format === :dense
        @test [f.value_type for f ∈ st.failed_attempts] == [Int64]
        @test all(f -> f.unit_pivots >= 0 && f.seconds >= 0, st.failed_attempts)
        @test sprint(show, st) == "ReductionStats(value_type = Int128, unit_pivots = 1, core = 1×1, failed = [Int64])"
        @test occursin("failed         Int64 after", sprint(show, MIME"text/plain"(), st))
        # the last type does not suffice: the attempt is recorded, then the error raised
        err = caught(() -> sed(B; value_type = Int64, stats = st))
        @test err isa CoefficientOverflowError
        @test st.value_type === nothing && [f.value_type for f ∈ st.failed_attempts] == [Int64]
        # an input entry that fits no type: no attempt, the statistics are reset all the same
        err = caught(() -> sed(B; value_type = Int32, stats = st))
        @test err isa CoefficientOverflowError && err.entry !== nothing
        @test st.value_type === nothing && isempty(st.failed_attempts)
        # reused: everything reset at each call
        @test sed(sparse([1 0; 0 1]); stats = st) == BigInt[1, 1]
        @test sprint(show, st) == "ReductionStats(value_type = Int32, unit_pivots = 2, core = 0×0)"
    end

    @testset "exported and public names" begin
        SED = SparseElementaryDivisors
        for name ∈ (:sparse_elementary_divisors, :ReductionConfig, :ReductionStats, :CoefficientOverflowError,
                    :read_matrix_market)
            @test Base.isexported(SED, name)
        end
        @static if VERSION >= v"1.11.0-DEV.469"
            for name ∈ (:DEFAULT_VALUE_TYPES, :ReductionTrace, :RowsData)
                @test Base.ispublic(SED, name) && !Base.isexported(SED, name)
            end
        end
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
        # a command that does not exist (the task that writes the script on the standard input
        # warns of the broken pipe: silenced, the logger is inherited by that task)
        Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            @test_throws ErrorException timed_magma_elementary_divisors(sparse([2 0; 0 3]);
                                                                        magma_command = "sed-no-such-magma-command")
        end
        # a fake Magma (a Julia script with a fixed answer, which checks the matrix it receives
        # on its standard input), for the whole path: script sent, output read and parsed
        exe = joinpath(Sys.BINDIR, Base.julia_exename())
        script = joinpath(mktempdir(), "fake_magma.jl")
        write(script, """
            s = read(stdin, String)
            if occursin("SparseMatrix(2, 2, [ 1, 1, 2, 1, 2, 3 ])", s) && occursin("ElementaryDivisors", s)
                print("Time: 0.010\\n##DIVS_BEGIN##\\n[ 1, 6 ]\\n##DIVS_END##\\n")
            else
                print("##MAGMA_ERROR_BEGIN##unexpected input##MAGMA_ERROR_END##")
            end
            """)
        if occursin(' ', exe) || occursin(' ', script)        # magma_command is split at spaces
            @info "Fake Magma test skipped: a space in the path of Julia or of the script."
        else
            cmd = "$exe --startup-file=no $script"
            A = [2 0; 0 3]
            mtx = joinpath(dirname(script), "A.mtx")
            write_matrix_market(mtx, A)
            for input ∈ ((sparse(A),), (A,), (Hecke.matrix(Hecke.ZZ, A),), (mtx,),
                         ([[(1, 2)], [(2, 3)]], 2, 2), (2, 2, [[(1, 2)], [(2, 3)]]))
                @test timed_magma_elementary_divisors(input...; magma_command = cmd) == (BigInt[1, 6], 0.010)
            end
            out = mktemp() do path, io
                redirect_stdout(io) do
                    @test magma_elementary_divisors(sparse(A); magma_command = cmd, verbose = true) == BigInt[1, 6]
                end
                flush(io)
                read(path, String)
            end
            @test occursin("Magma: ElementaryDivisors computed in 0.01 s", out)
        end
    end
end

include("test_fingerprints.jl")