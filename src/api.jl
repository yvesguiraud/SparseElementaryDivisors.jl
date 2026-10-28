# api.jl — public interface.

"""
    SparseElementaryDivisors.DEFAULT_VALUE_TYPES

The sequence of integer types tried by `value_type = :auto`, in increasing width:
`(Int32, Int64, Int128, BigInt)`. Int32 comes first because it saves time and memory on
large matrices (B*(G₃₇) in trivial coefficients: 20 % of the time, 17 % of the peak memory
against Int64), and its failed attempts cost little when the coefficients grow (their
overflow comes early).
"""
const DEFAULT_VALUE_TYPES = (Int32, Int64, Int128, BigInt)

# Signed integers of fixed width: those of Base, or wider ones such as Int256 of
# BitIntegers.jl (their *_with_overflow methods give the checked arithmetic).
_is_fixed_width(T) = T <: Signed && isbitstype(T)

# The sequence of types given by the keyword value_type, checked.
function _value_types(value_type)
    types = value_type === :auto ? DEFAULT_VALUE_TYPES :
            value_type isa Type ? (value_type,) :
            value_type isa Union{Tuple, AbstractVector} ? Tuple(value_type) :
            throw(ArgumentError("value_type must be :auto, an integer type or a sequence of " *
                                "integer types (got $(repr(value_type)))"))
    isempty(types) && throw(ArgumentError("value_type: empty sequence of types"))
    for (i, T) ∈ enumerate(types)
        (T isa Type && (T === BigInt || _is_fixed_width(T))) ||
            throw(ArgumentError("value_type: $(repr(T)) is neither a signed integer type of " *
                                "fixed width (Int8 to Int128, or wider ones such as Int256 of " *
                                "BitIntegers.jl) nor BigInt"))
        if i > 1
            P = types[i - 1]
            (P !== BigInt && (T === BigInt || sizeof(T) > sizeof(P))) ||
                throw(ArgumentError("value_type: the types must be of strictly increasing width, " *
                                    "BigInt last (got $P before $T)"))
        end
    end
    return types
end

# Does every input entry fit in T? `ext`: extreme entries returned by _check_rows.
function _fits(::Type{T}, ext) where {T}
    (T === BigInt || ext === nothing) && return true
    return typemin(T) <= ext.lo[3] && ext.hi[3] <= typemax(T)
end

# The error for an input entry that does not fit in T (one of the two extreme entries).
function _input_overflow(::Type{T}, ext) where {T}
    r, c, v = ext.hi[3] > typemax(T) ? ext.hi : ext.lo
    return CoefficientOverflowError(T, 0, (r, c, BigInt(v)))
end

# Check of the input given row by row, in one pass over the entries, before any
# construction: `DualSparseMatrix` writes with `@inbounds` and its merges assume strictly
# increasing column indices in each row. Stops at the first fault, with an ArgumentError
# in terms of the matrix (row, column), whatever the form of the input. Returns the
# extreme nonzero entries, `(lo = (row, column, value), hi = (row, column, value))`, used
# for the first type of the widening, or `nothing` if there is no nonzero entry.
function _check_rows(rows::RowsData{V}, m::Int, n::Int) where {V<:Integer}
    (m >= 0 && n >= 0) ||
        throw(ArgumentError("m and n must be nonnegative (got m = $m, n = $n)"))
    length(rows) == m ||
        throw(ArgumentError("rows has $(length(rows)) rows but m = $m"))
    found = false
    lo_r = lo_c = hi_r = hi_c = 0
    local lo_v::V, hi_v::V   # set at the first nonzero entry (V may be abstract: no zero(V))
    for r ∈ 1:m
        row = rows[r]
        prev = 0
        for k ∈ eachindex(row)
            c, v = row[k]
            1 <= c <= n ||
                throw(ArgumentError("row $r (entry $k): column $c is outside 1:$n"))
            if c <= prev
                what = c == prev ? "column $c appears twice" : "column $c after column $prev"
                throw(ArgumentError("row $r (entry $k): $what; the columns of each row must be " *
                                    "strictly increasing (sparse(I, J, V) sorts and sums duplicates)"))
            end
            prev = c
            iszero(v) && continue
            if !found || v < lo_v
                lo_v = v; lo_r = r; lo_c = c
            end
            if !found || v > hi_v
                hi_v = v; hi_r = r; hi_c = c
            end
            found = true
        end
    end
    return found ? (lo = (lo_r, lo_c, lo_v), hi = (hi_r, hi_c, hi_v)) : nothing
end

# Values of the rows (option A of s3): Int when every entry fits a machine word, BigInt
# otherwise (the conversions below then pass the type to a function barrier).
_fits_int(v) = typemin(Int) <= v <= typemax(Int)

# Conversion from CSC: rows as (column, value) pairs sorted by increasing column.
_rows_data(A::SparseMatrixCSC{<:Integer}) = _rows_data(A, all(_fits_int, nonzeros(A)) ? Int : BigInt)

function _rows_data(A::SparseMatrixCSC{<:Integer}, ::Type{V}) where {V}
    m, n = size(A)
    rows = [Tuple{Int,V}[] for _ ∈ 1:m]
    rv = rowvals(A)
    nz = nonzeros(A)
    for j ∈ 1:n, k ∈ nzrange(A, j)
        v = nz[k]
        iszero(v) && continue
        push!(rows[rv[k]], (j, V(v)))
    end
    return rows
end

# Conversion from a Hecke sparse matrix.
_rows_data(S::Hecke.SMat{Hecke.ZZRingElem}) =
    _rows_data(S, all(row -> all(_fits_int, row.values), S.rows) ? Int : BigInt)

function _rows_data(S::Hecke.SMat{Hecke.ZZRingElem}, ::Type{V}) where {V}
    rows = Vector{Vector{Tuple{Int,V}}}(undef, Hecke.nrows(S))
    for (i, row) ∈ enumerate(S.rows)
        r = Tuple{Int,V}[]
        sizehint!(r, length(row.pos))
        for (c, v) ∈ zip(row.pos, row.values)
            iszero(v) || push!(r, (Int(c), V(v)))
        end
        rows[i] = r
    end
    return rows
end

# Conversion from a Flint dense matrix.
function _rows_data(A::Hecke.ZZMatrix)
    m, n = Hecke.nrows(A), Hecke.ncols(A)
    V = all(_fits_int(A[i, j]) for i ∈ 1:m, j ∈ 1:n) ? Int : BigInt
    return _rows_data(A, V)
end

function _rows_data(A::Hecke.ZZMatrix, ::Type{V}) where {V}
    m, n = Hecke.nrows(A), Hecke.ncols(A)
    rows = Vector{Vector{Tuple{Int,V}}}(undef, m)
    for i ∈ 1:m
        r = Tuple{Int,V}[]
        for j ∈ 1:n
            v = A[i, j]
            iszero(v) || push!(r, (j, V(v)))
        end
        rows[i] = r
    end
    return rows
end

"""
    sparse_elementary_divisors(A; value_type = :auto, config = ReductionConfig(), stats = nothing,
                               verbose = false, settings...)
    sparse_elementary_divisors(rows, m, n; kwargs...)
    sparse_elementary_divisors(m, n, rows; kwargs...)
    sparse_elementary_divisors(path::AbstractString; kwargs...)

Non-zero elementary divisors of the integer matrix `A`, in increasing order: a `Vector{BigInt}` with `d[i] ∣ d[i+1]`.

`A` can be a `SparseMatrixCSC` (fast path), any other `AbstractSparseMatrix` or `AbstractMatrix{<:Integer}` (converted), a Hecke matrix (`Hecke.SMat{ZZRingElem}` or `Hecke.ZZMatrix`), or the matrix given row by row: `rows[r]` is the vector of `(column, value)` pairs of row `r`, with columns in `1:n` and strictly increasing, for an `m × n` matrix (zero values are allowed and ignored; values of any integer type, entries beyond a machine word included). The input is checked before the computation: an `ArgumentError` names the first fault (column out of range, unsorted or repeated column, wrong number of rows); unsorted triplets go through `sparse(I, J, V)`, which sorts them and sums duplicates. With a file name, the matrix is read in Matrix Market format (see [`read_matrix_market`](@ref)).

- `value_type`: integer type(s) of the coefficients during the elimination, with checked arithmetic. `:auto` (default) tries the types of `SparseElementaryDivisors.DEFAULT_VALUE_TYPES` = `(Int32, Int64, Int128, BigInt)` in turn, starting with the first one that holds all the entries, and restarts the computation from the input with the next type after an overflow (each matrix independently: a failed attempt costs its time). A type (`value_type = Int64`) uses that type only; a sequence (`value_type = (Int64, Int128)`) tries its types in turn. Accepted types: signed integers of fixed width (`Int8` to `Int128`, or wider ones such as `Int256` of BitIntegers.jl) in strictly increasing width, and `BigInt` last. When the last type does not suffice, a [`CoefficientOverflowError`](@ref) is raised.
- `config`: the settings of the reduction, a [`ReductionConfig`](@ref) (pivot search, parallelism, format of the core, memory); the same for all the attempts of the widening.
- `settings...`: any field of `ReductionConfig` given directly as a keyword, which takes precedence over `config`. The most common: `nworkers` (number of parallel tasks, default `Threads.nthreads()`; start Julia with `-t auto` or `-t N`; the result does not depend on it) and `core_format` (`:auto`, `:sparse` or `:dense`, format of the core passed to Hecke). An unknown keyword raises an `ArgumentError`.
- `stats`: a [`ReductionStats`](@ref), reset then filled in place: type retained, failed attempts, counters, times and core.
- `verbose`: print the attempts that overflowed, the type retained, and the size and format of the core.

Developer tools: `trace`, a `SparseElementaryDivisors.ReductionTrace` (trajectory fingerprint of the attempt retained), and `probe`, a function of no argument called after each Markowitz pivot.
"""
function sparse_elementary_divisors(rows::RowsData{<:Integer}, m::Int, n::Int;
                                    value_type = :auto,
                                    config::ReductionConfig = ReductionConfig(),
                                    stats::Union{Nothing, ReductionStats} = nothing,
                                    trace::Union{Nothing, ReductionTrace} = nothing,
                                    probe = nothing,
                                    verbose::Bool = false,
                                    settings...)
    stats === nothing || _reset!(stats)
    types = _value_types(value_type)
    cfg = _with_settings(config, settings)
    ext = _check_rows(rows, m, n)
    first_type = findfirst(T -> _fits(T, ext), types)
    first_type === nothing && throw(_input_overflow(last(types), ext))
    d = _reduce_and_divisors_widening(rows, m, n, types[first_type:end]; config = cfg, trace, stats,
                                      probe, verbose)
    return BigInt.(d)
end

# Form (m, n, rows), with the same order of objects as in the saved test files.
sparse_elementary_divisors(m::Int, n::Int, rows::RowsData{<:Integer}; kwargs...) =
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
