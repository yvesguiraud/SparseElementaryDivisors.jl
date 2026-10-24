# Trajectory fingerprints: non-regression test of the reduction.
#
# The references in data/references.toml were computed with SparseElementaryDivisors
# v0.1.0 (see data/README.md). For every matrix, the elementary divisors and the
# trajectory fingerprint (ReductionTrace) must be identical with the default
# configuration and with the parallel code paths forced (search_mode = :always,
# apply_min_work = 0), for 1, 2 and 4 workers. Matrices that overflowed with a smaller
# value type are also checked to raise the overflow error with that type.
#
# Heavy matrices (G37) are skipped unless the environment variable SED_HEAVY_DIR gives
# the folder that contains their .mtx files.

using TOML

const FP_DATA = joinpath(@__DIR__, "data")
const FP_REFS = TOML.parsefile(joinpath(FP_DATA, "references.toml"))
const FP_HEAVY_DIR = get(ENV, "SED_HEAVY_DIR", "")
const FP_VALUE_TYPES = Dict("Int32" => Int32, "Int64" => Int64, "BigInt" => BigInt)

function fp_torsion(d)
    t = Dict{BigInt, Int}()
    for x in d
        x > 1 && (t[x] = get(t, x, 0) + 1)
    end
    return ["$k:$(t[k])" for k in sort!(collect(keys(t)))]
end

fp_hex(h) = "0x" * string(h; base = 16, pad = 16)

function fp_check(ref, path)
    SED = SparseElementaryDivisors
    A = read_matrix_market(path)
    m, n = size(A)
    @test (m, n) == (ref["rows"], ref["cols"])
    @test nnz(A) == ref["nnz"]
    rows = SED._rows_data(A)

    for t in ref["overflow"]
        err = try
            SED._reduce_and_divisors(rows, m, n; value_type = FP_VALUE_TYPES[t])
            nothing
        catch e
            e
        end
        @test err isa CoefficientOverflowError && err.value_type === FP_VALUE_TYPES[t]
    end

    vt = FP_VALUE_TYPES[ref["value_type"]]
    for (nw, forced) in ((1, false), (1, true), (2, true), (4, true))
        tr = SED.ReductionTrace()
        cfg = forced ?
            SED.ReductionConfig(; nworkers = nw, search_mode = :always, apply_min_work = 0, trace = tr) :
            SED.ReductionConfig(; nworkers = nw, trace = tr)
        d = SED._reduce_and_divisors(rows, m, n; value_type = vt, config = cfg)
        @test length(d) == ref["rank"]
        @test count(==(1), d) == ref["units"]
        @test fp_torsion(d) == ref["torsion"]
        @test fp_hex(tr.h) == ref["trace_hash"]
        @test (tr.n_pivots, tr.n_sge, tr.sum_ones) ==
              (ref["trace_pivots"], ref["trace_sge"], ref["trace_sge_units"])
    end
end

# The widening (step 3 of v0.1.1): a matrix that overflows Int32 gives, from the default
# sequence of types, the divisors and the trajectory fingerprint of its reference type.
function fp_check_widening(ref, path)
    SED = SparseElementaryDivisors
    A = read_matrix_market(path)
    m, n = size(A)
    rows = SED._rows_data(A)
    vt = FP_VALUE_TYPES[ref["value_type"]]
    for nw in (1, 4)
        tr = SED.ReductionTrace()
        st = SED.ReductionStats()
        cfg = SED.ReductionConfig(; nworkers = nw, trace = tr, stats = st)
        d = SED._reduce_and_divisors_widening(rows, m, n, SED.DEFAULT_VALUE_TYPES; config = cfg)
        @test st.value_type === vt
        @test [a.value_type for a in st.failed_attempts] == [FP_VALUE_TYPES[t] for t in ref["overflow"]]
        @test all(a -> a.unit_pivots >= 0 && a.seconds >= 0, st.failed_attempts)
        @test length(d) == ref["rank"]
        @test count(==(1), d) == ref["units"]
        @test fp_torsion(d) == ref["torsion"]
        @test fp_hex(tr.h) == ref["trace_hash"]
        @test (tr.n_pivots, tr.n_sge, tr.sum_ones) ==
              (ref["trace_pivots"], ref["trace_sge"], ref["trace_sge_units"])
    end
    # the public entry point gives the same divisors without any setting
    d = sparse_elementary_divisors(A)
    @test length(d) == ref["rank"] && count(==(1), d) == ref["units"] && fp_torsion(d) == ref["torsion"]
end

@testset "trajectory fingerprints" begin
    for ref in FP_REFS["matrix"]
        heavy = ref["heavy"]
        heavy && isempty(FP_HEAVY_DIR) && continue
        dir = heavy ? FP_HEAVY_DIR : FP_DATA
        @testset "$(ref["file"])" begin
            fp_check(ref, joinpath(dir, ref["file"]))
        end
    end
    isempty(FP_HEAVY_DIR) &&
        @info "Heavy fingerprint tests (G37) skipped; set SED_HEAVY_DIR to run them."
    for ref in FP_REFS["matrix"]
        (isempty(ref["overflow"]) || ref["heavy"]) && continue
        @testset "widening $(ref["file"])" begin
            fp_check_widening(ref, joinpath(FP_DATA, ref["file"]))
        end
    end
end
