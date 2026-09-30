# api.jl — public interface.

const RowsData = Vector{Vector{Tuple{Int,Int}}}

function _check_int(v, i, j)
    typemin(Int) <= v <= typemax(Int) ||
        throw(ArgumentError("entry ($i, $j) = $v is outside the range of Int (BigInt entries are not supported)"))
    return Int(v)
end

# Conversion from CSC: rows as (column, value) pairs sorted by increasing column.
function _rows_data(A::SparseMatrixCSC{<:Integer})
    m, n = size(A)
    rows = [Tuple{Int,Int}[] for _ ∈ 1:m]
    rv = rowvals(A)
    nz = nonzeros(A)
    for j ∈ 1:n, k ∈ nzrange(A, j)
        v = nz[k]
        iszero(v) && continue
        push!(rows[rv[k]], (j, _check_int(v, rv[k], j)))
    end
    return rows
end

# Conversion from a Hecke sparse matrix.
function _rows_data(S::Hecke.SMat{Hecke.ZZRingElem})
    rows = Vector{Vector{Tuple{Int,Int}}}(undef, Hecke.nrows(S))
    for (i, row) ∈ enumerate(S.rows)
        r = Tuple{Int,Int}[]
        sizehint!(r, length(row.pos))
        for (c, v) ∈ zip(row.pos, row.values)
            iszero(v) || push!(r, (Int(c), _check_int(v, i, c)))
        end
        rows[i] = r
    end
    return rows
end

# Conversion from a Flint dense matrix.
function _rows_data(A::Hecke.ZZMatrix)
    m, n = Hecke.nrows(A), Hecke.ncols(A)
    rows = Vector{Vector{Tuple{Int,Int}}}(undef, m)
    for i ∈ 1:m
        r = Tuple{Int,Int}[]
        for j ∈ 1:n
            v = A[i, j]
            iszero(v) || push!(r, (j, _check_int(v, i, j)))
        end
        rows[i] = r
    end
    return rows
end

"""
    sparse_elementary_divisors(A; value_type = Int32, nworkers = nothing, core_format = :auto,
                               acceptable_score_threshold = 4, max_buckets_scanned = 1,
                               verbose = false)
    sparse_elementary_divisors(rows, m, n; kwargs...)
    sparse_elementary_divisors(m, n, rows; kwargs...)
    sparse_elementary_divisors(path::AbstractString; kwargs...)

Non-zero elementary divisors of the integer matrix `A`, in increasing order: a `Vector{BigInt}` with `d[i] ∣ d[i+1]`.

`A` can be a `SparseMatrixCSC` (fast path), any other `AbstractSparseMatrix` or `AbstractMatrix{<:Integer}` (converted), a Hecke matrix (`Hecke.SMat{ZZRingElem}` or `Hecke.ZZMatrix`), or the matrix given row by row: `rows[r]` is the list of `(column, value)` pairs of row `r`, sorted by increasing column, for an `m × n` matrix. With a file name, the matrix is read in Matrix Market format (see [`read_matrix_market`](@ref)).

- `value_type`: integer type used during the elimination (`Int32` by default, the fastest). An arithmetic overflow raises an explicit `ErrorException`; retry then with `Int64`, `Int128` or `BigInt`.
- `nworkers`: number of parallel tasks (default: `Threads.nthreads()`); start Julia with `-t auto` or `-t N` to use them. The result does not depend on `nworkers`.
- `core_format`: `:auto`, `:sparse` or `:dense`, format of the core passed to Hecke.
- `acceptable_score_threshold`, `max_buckets_scanned`: parameters of the pivot search.
- `verbose`: print the size and format of the core.
"""
function sparse_elementary_divisors(rows::RowsData, m::Int, n::Int;
                                    value_type::Type = Int32,
                                    nworkers::Union{Nothing,Integer} = nothing,
                                    core_format::Symbol = :auto,
                                    acceptable_score_threshold::Int = 4,
                                    max_buckets_scanned::Int = 1,
                                    verbose::Bool = false)
    # The configuration is built at run time (not at precompilation time): the default `nworkers`
    # depends on `Threads.nthreads()` of the current session.
    config = nworkers === nothing ? ReductionConfig() : ReductionConfig(; nworkers = Int(nworkers))
    d = _reduce_and_divisors(rows, m, n; value_type, config, core_format,
                             acceptable_score_threshold, max_buckets_scanned, verbose)
    return BigInt.(d)
end

# Form (m, n, rows), with the same order of objects as in the saved test files.
sparse_elementary_divisors(m::Int, n::Int, rows::RowsData; kwargs...) =
    sparse_elementary_divisors(rows, m, n; kwargs...)

function sparse_elementary_divisors(A::SparseMatrixCSC{<:Integer}; kwargs...)
    m, n = size(A)
    return sparse_elementary_divisors(_rows_data(A), m, n; kwargs...)
end

# Other sparse matrix types (from other packages): converted to SparseMatrixCSC.
sparse_elementary_divisors(A::AbstractSparseMatrix{<:Integer}; kwargs...) =
    sparse_elementary_divisors(SparseMatrixCSC(A); kwargs...)

# Dense matrices and wrappers (`Adjoint`, `Transpose`, views…): converted to sparse.
sparse_elementary_divisors(A::AbstractMatrix{<:Integer}; kwargs...) =
    sparse_elementary_divisors(sparse(A); kwargs...)

# Hecke matrices.
function sparse_elementary_divisors(S::Union{Hecke.SMat{Hecke.ZZRingElem},Hecke.ZZMatrix}; kwargs...)
    return sparse_elementary_divisors(_rows_data(S), Hecke.nrows(S), Hecke.ncols(S); kwargs...)
end

# Matrix Market file.
sparse_elementary_divisors(path::AbstractString; kwargs...) =
    sparse_elementary_divisors(read_matrix_market(path); kwargs...)
