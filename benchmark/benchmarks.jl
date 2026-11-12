# Benchmarks of SparseElementaryDivisors.jl, outside the package: times of the package, of Hecke alone and of Magma on the matrices of the page Benchmarks of the documentation.
#
# The matrices (G37_d3.mtx to G37_d6.mtx, ZT2w4m5_d2.mtx to ZT2w4m5_d4.mtx) are not in the repository: the environment variable SED_BENCHMARK_DIR gives their folder. One process per tool and matrix, from the root of the repository:
#
#     julia --project=benchmark -e "using Pkg; Pkg.instantiate()"        (the first time)
#     julia -t auto --project=benchmark benchmark/benchmarks.jl sed G37_d4
#     julia -t auto --project=benchmark benchmark/benchmarks.jl sed ZT2w4m5_d4 runs=2
#     julia --project=benchmark benchmark/benchmarks.jl hecke G37_d4 runs=1
#     julia --project=benchmark benchmark/benchmarks.jl magma G37_d4
#
# Tools:
# - `sed`, the package: after a warm-up on G34_d3 of test/data, k runs with one worker and k runs with Threads.nthreads() workers, interleaved (option workers = 1, all or both, the default);
# - `hecke`, Hecke.elementary_divisors on the matrix as a Hecke.SMat, built beforehand;
# - `magma`, Magma's ElementaryDivisors on this machine, with the time measured by Magma for that call alone.
# The reading of the file and the construction of the SMat are not timed; a full garbage collection precedes each run. Option runs = k, 3 by default.
#
# Each run adds a line to benchmark/results/<processor>.csv: machine, versions, tool, workers, matrix, time, statistics of the reduction (sed), peak memory of the process so far, and a summary of the divisors. The divisors must be identical from one run to the next, and equal to the references of test/data/references.toml when the matrix has one (G37).

using SparseElementaryDivisors, SparseArrays, TOML
import Hecke

const SED = SparseElementaryDivisors
const ROOT = pkgdir(SED)
const RESULTS = joinpath(@__DIR__, "results")

const HEADER = ["date", "cpu", "memory_gib", "os", "julia_threads", "julia", "sed", "hecke", "tool", "nworkers", "matrix", "rows", "cols", "nnz", "run", "seconds", "value_type", "failed", "unit_pivots", "core_rows", "core_cols", "core_nnz", "core_format", "core_seconds", "maxrss_gib", "divisors"]

# ---------------------------------------------------------------------------
# Machine and versions

cpu_model() = strip(replace(Sys.cpu_info()[1].model, r"\s+" => " "))

# A file name from the processor: "Intel(R) Core(TM) Ultra 9 285" gives "intel-core-ultra-9-285".
function machine_slug(model::AbstractString)
    s = replace(lowercase(model), "(r)" => "", "(tm)" => "", r"\bcpu\b" => "", r"@.*$" => "")
    return strip(replace(strip(s), r"[^a-z0-9]+" => "-"), '-')
end

os_name() = Sys.iswindows() ? "Windows" : Sys.islinux() ? "Linux" : Sys.isapple() ? "macOS" : string(Sys.KERNEL)

# Version of the package, with the commit of the repository ("+" when tracked files other than the results differ from it).
function sed_version()
    commit = try
        c = readchomp(`git -C $ROOT rev-parse --short HEAD`)
        changes = readchomp(`git -C $ROOT status --porcelain --untracked-files=no -- . ":(exclude)benchmark/results"`)
        isempty(changes) ? c : c * "+"
    catch
        "unknown"
    end
    return string(pkgversion(SED), " ", commit)
end

function machine_fields()
    return (date = Libc.strftime("%Y-%m-%d %H:%M", time()), cpu = cpu_model(), memory_gib = round(Sys.total_memory() / 2^30; digits = 1), os = os_name(), julia_threads = Threads.nthreads(), julia = string(VERSION), sed = sed_version(), hecke = string(pkgversion(Hecke)))
end

# ---------------------------------------------------------------------------
# Results

csv_field(x) = (s = string(x); occursin(r"[\",\n]", s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s)

function append_row(path::AbstractString, row::NamedTuple)
    keys(row) == Tuple(Symbol.(HEADER)) || error("internal: the fields of the row do not match the header")
    newfile = !isfile(path)
    open(path, "a") do io
        newfile && println(io, join(HEADER, ","))
        println(io, join(map(csv_field, values(row)), ","))
    end
    return nothing
end

# "1^41221 2 6": the number of divisors equal to 1, then the others.
divisors_summary(d) = join(["1^$(count(isone, d))"; [string(x) for x in d if !isone(x)]], " ")

# The divisors greater than 1 with their multiplicities, in the form of test/data/references.toml.
function torsion_strings(d)
    t = Dict{BigInt,Int}()
    for x in d
        x > 1 && (t[x] = get(t, x, 0) + 1)
    end
    return ["$k:$(t[k])" for k in sort!(collect(keys(t)))]
end

function check_reference(name::AbstractString, d)
    refs = TOML.parsefile(joinpath(ROOT, "test", "data", "references.toml"))["matrix"]
    i = findfirst(r -> r["file"] == name * ".mtx", refs)
    i === nothing && return nothing
    ref = refs[i]
    length(d) == ref["rank"] && count(isone, d) == ref["units"] && torsion_strings(d) == ref["torsion"] ||
        error("$name: the divisors differ from those of test/data/references.toml")
    return nothing
end

# The divisors of the first run, against which the following ones are checked.
function check_same!(first_divisors::Vector{Vector{BigInt}}, name::AbstractString, d)
    if isempty(first_divisors)
        check_reference(name, d)
        push!(first_divisors, BigInt.(d))
    elseif d != first_divisors[1]
        error("$name: the divisors differ from those of the first run")
    end
    return nothing
end

function base_row(machine, tool, nworkers, name, A, run, seconds, d)
    m, n = size(A)
    return (; machine..., tool, nworkers, matrix = name, rows = m, cols = n, nnz = nnz(A), run, seconds = round(seconds; digits = 3), value_type = "", failed = "", unit_pivots = "", core_rows = "", core_cols = "", core_nnz = "", core_format = "", core_seconds = "", maxrss_gib = round(Sys.maxrss() / 2^30; digits = 2), divisors = divisors_summary(d))
end

function report(tool, nworkers, name, run, seconds)
    w = nworkers == "" ? "" : ", $nworkers workers"
    println("$tool $name$w, run $run: $(round(seconds; digits = 2)) s")
    flush(stdout)
end

# ---------------------------------------------------------------------------
# Tools

function run_sed(machine, csv, name, A, runs, workers)
    N = Threads.nthreads()
    nws = workers == "1" ? [1] : workers == "all" ? [N] : workers == "both" ? unique([1, N]) : error("workers = 1, all or both")
    W = read_matrix_market(joinpath(ROOT, "test", "data", "G34_d3.mtx"))
    for nw in nws
        sparse_elementary_divisors(W; nworkers = nw)
    end
    first_divisors = Vector{BigInt}[]
    for run in 1:runs, nw in nws
        st = ReductionStats()
        GC.gc()
        t0 = time_ns()
        d = sparse_elementary_divisors(A; nworkers = nw, stats = st)
        seconds = (time_ns() - t0) / 1e9
        check_same!(first_divisors, name, d)
        row = merge(base_row(machine, "sed", nw, name, A, run, seconds, d), (value_type = string(st.value_type), failed = join((string(a.value_type) for a in st.failed_attempts), " "), unit_pivots = st.unit_pivots, core_rows = st.core_size.rows, core_cols = st.core_size.cols, core_nnz = st.core_size.nnz, core_format = string(st.core_format), core_seconds = round(st.core_seconds; digits = 3)))
        append_row(csv, row)
        report("sed", nw, name, run, seconds)
    end
    return nothing
end

# The matrix as a Hecke.SMat, row by row (as in the earlier measurements of Hecke alone).
function hecke_matrix(A::SparseMatrixCSC)
    m, n = size(A)
    At = copy(transpose(A))   # column i of At is row i of A
    rv = rowvals(At)
    nz = nonzeros(At)
    S = Hecke.sparse_matrix(Hecke.ZZ, 0, n)
    for i in 1:m
        entries = Tuple{Int,Int}[(rv[p], Int(nz[p])) for p in nzrange(At, i) if !iszero(nz[p])]
        push!(S, Hecke.sparse_row(Hecke.ZZ, entries))
    end
    return S
end

function run_hecke(machine, csv, name, A, runs)
    first_divisors = Vector{BigInt}[]
    for run in 1:runs
        S = hecke_matrix(A)
        GC.gc()
        t0 = time_ns()
        e = Hecke.elementary_divisors(S)
        seconds = (time_ns() - t0) / 1e9
        d = sort!([BigInt(x) for x in e if !iszero(x)])
        check_same!(first_divisors, name, d)
        append_row(csv, base_row(machine, "hecke", "", name, A, run, seconds, d))
        report("hecke", "", name, run, seconds)
    end
    return nothing
end

function run_magma(machine, csv, name, A, runs)
    first_divisors = Vector{BigInt}[]
    for run in 1:runs
        GC.gc()
        d, seconds = timed_magma_elementary_divisors(A)
        check_same!(first_divisors, name, d)
        append_row(csv, base_row(machine, "magma", "", name, A, run, seconds, d))
        report("magma", "", name, run, seconds)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Command line

function main(args)
    length(args) >= 2 || error("usage: benchmarks.jl sed|hecke|magma <matrix> [runs=k] [workers=1|all|both]")
    tool = args[1]
    name = args[2]
    opts = Dict{String,String}()
    for a in args[3:end]
        occursin('=', a) || error("option $a: expected key=value")
        k, v = split(a, '='; limit = 2)
        opts[k] = v
    end
    runs = parse(Int, get(opts, "runs", "3"))
    workers = get(opts, "workers", "both")
    dir = get(ENV, "SED_BENCHMARK_DIR", "")
    isempty(dir) && error("SED_BENCHMARK_DIR must give the folder of the .mtx files")
    A = read_matrix_market(joinpath(dir, name * ".mtx"))
    machine = machine_fields()
    mkpath(RESULTS)
    csv = joinpath(RESULTS, machine_slug(machine.cpu) * ".csv")
    println("$(machine.cpu), $(machine.julia_threads) threads; SparseElementaryDivisors $(machine.sed), Hecke $(machine.hecke); results in $csv")
    if tool == "sed"
        run_sed(machine, csv, name, A, runs, workers)
    elseif tool == "hecke"
        run_hecke(machine, csv, name, A, runs)
    elseif tool == "magma"
        run_magma(machine, csv, name, A, runs)
    else
        error("tool $tool: expected sed, hecke or magma")
    end
    return nothing
end

main(ARGS)
