# solver.jl — core reduction: sparse elimination with unit (±1) pivots, Markowitz-style pivot choice,
# parallel pivot search and fill-in application, exact checked integer arithmetic.

"""
    RowsData{V}

A matrix given row by row: `rows[r]` is the vector of the `(column, value)` pairs of row `r`,
columns strictly increasing, values of type `V <: Integer` (`Int`, or `BigInt` for entries
beyond a machine word).
"""
const RowsData{V} = Vector{Vector{Tuple{Int, V}}}

######################
### Degree buckets ###
######################

mutable struct DegreeBuckets
    buckets::Vector{Vector{Int}}   # buckets[d+1] = active indices of degree d, for d ∈ 0:max_deg
                                   # buckets[max_deg+2] = overflow bucket, degree > max_deg
    pos_in_bucket::Vector{Int}     # current position of each index in its current bucket
    bucket_of::Vector{Int}         # (1-based) index of the current bucket of each entity
    max_deg::Int
    low::Int                       # cursor (P5): the buckets of the degrees d with 2 ≤ d < low are
                                   # empty; lowered by _bucket_insert!, raised by the pivot search
end

function DegreeBuckets(n::Int, initial_degrees::Vector{Int}, max_deg::Int)
    buckets = [Int[] for _ ∈ 1:(max_deg + 2)]
    pos_in_bucket = zeros(Int, n)
    bucket_of = zeros(Int, n)
    db = DegreeBuckets(buckets, pos_in_bucket, bucket_of, max_deg, 2)
    for i ∈ 1:n
        bidx = degree_to_bucket_idx(db, initial_degrees[i])
        bucket = db.buckets[bidx]
        push!(bucket, i)
        db.pos_in_bucket[i] = length(bucket)
        db.bucket_of[i] = bidx
    end
    return db
end

@inline degree_to_bucket_idx(db::DegreeBuckets, d::Int)::Int =
    d <= db.max_deg ? d + 1 : db.max_deg + 2

# The only way, besides the constructor (where `low = 2` holds trivially), for an index to
# enter a bucket: keeps the invariant of `low`.
function _bucket_insert!(db::DegreeBuckets, i::Int, d::Int)
    bidx = degree_to_bucket_idx(db, d)
    bucket = db.buckets[bidx]
    push!(bucket, i)
    db.pos_in_bucket[i] = length(bucket)
    db.bucket_of[i] = bidx
    2 <= d < db.low && (db.low = d)
    return nothing
end

function _bucket_remove!(db::DegreeBuckets, i::Int)
    bidx = db.bucket_of[i]
    bidx == 0 && return nothing
    bucket = db.buckets[bidx]
    p = db.pos_in_bucket[i]
    last_elem = bucket[end]
    bucket[p] = last_elem
    db.pos_in_bucket[last_elem] = p
    pop!(bucket)
    db.bucket_of[i] = 0
    db.pos_in_bucket[i] = 0
    return nothing
end

"""
    update_degrees!(db, i, new_d)

Move entity `i` to the bucket matching its new degree `new_d`, if the degree changed bucket. Amortized O(1) cost.
"""
function update_degrees!(db::DegreeBuckets, i::Int, new_d::Int)
    new_bidx = degree_to_bucket_idx(db, new_d)
    db.bucket_of[i] == new_bidx && return nothing
    _bucket_remove!(db, i)
    _bucket_insert!(db, i, new_d)
    return nothing
end

###########################################
### Index of unit values (generic in T) ###
###########################################

"""
    is_unit(v::T) where T

True if `v` equals 1 or -1 in type `T`.
"""
@inline is_unit(v::T) where {T} = v == one(T) || v == -one(T)

mutable struct UnitColumnIndex
    unit_count::Vector{Int}    # number of ±1 entries, per column
    buckets::DegreeBuckets      # subset of col_buckets restricted to unit_count[c] > 0
end

"""
    UnitColumnIndex(col_rows, col_vals, n, max_deg)

Initial construction by a full scan (done once, when the matrix is created): all later updates are incremental (see remove_entry! and merge_col!).
"""
function UnitColumnIndex(col_rows::Vector{Vector{Int}}, col_vals::Vector{Vector{T}},
                          n::Int, max_deg::Int) where {T}
    unit_count = zeros(Int, n)
    buckets = [Int[] for _ ∈ 1:(max_deg + 2)]
    pos_in_bucket = zeros(Int, n)
    bucket_of = zeros(Int, n)
    db = DegreeBuckets(buckets, pos_in_bucket, bucket_of, max_deg, 2)
    for c ∈ 1:n
        cnt = 0
        @inbounds for v ∈ col_vals[c]
            is_unit(v) && (cnt += 1)
        end
        unit_count[c] = cnt
        cnt > 0 && _bucket_insert!(db, c, length(col_rows[c]))
    end
    return UnitColumnIndex(unit_count, db)
end

function unit_index_sync!(ui::UnitColumnIndex, c::Int, new_count::Int, new_deg::Int)
    old_count = ui.unit_count[c]
    ui.unit_count[c] = new_count
    if old_count == 0 && new_count > 0
        _bucket_insert!(ui.buckets, c, new_deg)
    elseif old_count > 0 && new_count == 0
        _bucket_remove!(ui.buckets, c)
    elseif old_count > 0 && new_count > 0
        update_degrees!(ui.buckets, c, new_deg)
    end
    return nothing
end

function unit_index_clear_column!(ui::UnitColumnIndex, c::Int)
    if ui.unit_count[c] > 0
        _bucket_remove!(ui.buckets, c)
        ui.unit_count[c] = 0
    end
    return nothing
end

########################################################################
### Generic matrix structure DualSparseMatrix{T} and its constructor ###
########################################################################

# Values (`row_vals`, `col_vals`) and value buffers are stored as `Vector{T}` / `Vector{Vector{T}}`; the merge output buffers `_merge_cols_buf` / `_merge_vals_buf` are one scratch buffer shared by row and column merges.

struct DualSparseMatrix{T<:Integer}
    row_cols::Vector{Vector{Int}}     # Sorted list of active columns per row (indices, Int)
    row_vals::Vector{Vector{T}}       # Corresponding values per row (coefficients, T)
    col_rows::Vector{Vector{Int}}     # Sorted list of active rows per column (indices, Int)
    col_vals::Vector{Vector{T}}       # Corresponding values per column (coefficients, T)
    num_rows::Int
    num_cols::Int

    # Pre-allocated work buffers (avoid pressure on the garbage collector).
    _sister_rows_buf::Vector{Int}      # sister rows (col_rows[p_c] \ {p_r}), sorted — indices
    _sister_vals_buf::Vector{T}        # corresponding values of col_vals[p_c] — coefficients
    _p_row_cols_buf::Vector{Int}       # columns of the pivot row (row_cols[p_r] \ {p_c}), sorted
    _p_row_vals_buf::Vector{T}         # corresponding values — coefficients
    _factors_buf::Vector{T}            # pivot factor for each sister row — coefficient
    _merge_cols_buf::Vector{Int}       # reusable output buffer for merges (indices)
    _merge_vals_buf::Vector{T}         # (shared sequentially by the two merge phases)

    row_buckets::DegreeBuckets
    col_buckets::DegreeBuckets

    unit_index::UnitColumnIndex        # index of the columns with a ±1 entry
end

# Degrees above this bound share the overflow bucket.
const _MAX_TRACKED_DEGREE = 10_000

# Largest degree with its own bucket (P5): a row has at most n entries and a column at most m,
# so that with max(m, n) no entity reaches the overflow bucket, as with the bound 10 000 when
# max(m, n) ≤ 10 000 (the buckets beyond max(m, n) were always empty): same trajectory, bit for
# bit. At least 1, so that the bucket of degree 1 (structured elimination) is never the overflow
# bucket.
_max_tracked_degree(m::Int, n::Int) = clamp(max(m, n), 1, _MAX_TRACKED_DEGREE)

"""
    DualSparseMatrix(rows_data, m, n; max_tracked_degree = _max_tracked_degree(m, n),
                     value_type::Type{T} = Int32)

Build the matrix in linear time from the vector of `(j, v)` pairs of each row, already sorted by increasing row and column number. `value_type` sets the storage type of the coefficients (`T = Int32` by default). `max_tracked_degree`: largest degree with its own bucket, `clamp(max(m, n), 1, 10_000)` by default.

Precondition, not checked here (the loop writes with `@inbounds` and the merges assume it): `length(rows_data) == m`, in each row columns in `1:n` and strictly increasing, and values that fit in `T`. The public entry point checks the first two once with `_check_rows` and chooses `T` accordingly; the conversions `_rows_data` produce such rows.
"""
function DualSparseMatrix(rows_data::RowsData{<:Integer}, m::Int, n::Int;
                            max_tracked_degree::Int = _max_tracked_degree(m, n),
                            value_type::Type{T} = Int32) where {T}
    row_cols = [Int[] for _ ∈ 1:m]
    row_vals = [T[] for _ ∈ 1:m]
    col_rows = [Int[] for _ ∈ 1:n]
    col_vals = [T[] for _ ∈ 1:n]

    for r ∈ 1:m
        @inbounds for idx_t ∈ eachindex(rows_data[r])
            @inbounds t = rows_data[r][idx_t]
            c = t[1]
            v = t[2]
            v == 0 && continue
            vt = T(v)

            push!(row_cols[r], c)
            push!(row_vals[r], vt)

            push!(col_rows[c], r)
            push!(col_vals[c], vt)
        end
    end

    max_dim = max(m, n)
    row_buckets = DegreeBuckets(m, length.(row_cols), max_tracked_degree)
    col_buckets = DegreeBuckets(n, length.(col_rows), max_tracked_degree)
    unit_index = UnitColumnIndex(col_rows, col_vals, n, max_tracked_degree)
    return DualSparseMatrix{T}(
        row_cols, row_vals, col_rows, col_vals, m, n,
        sizehint!(Int[], max_dim), sizehint!(T[], max_dim),
        sizehint!(Int[], max_dim), sizehint!(T[], max_dim),
        T[],
        sizehint!(Int[], max_dim), sizehint!(T[], max_dim),
        row_buckets, col_buckets,
        unit_index
    )
end

########################################
### Point removal of an entry (SGE) ###
########################################

# Used only by the structured elimination (below), which removes entries as they are: no
# arithmetic, hence no overflow risk here. Removes the entry (r, c) if present, on the column
# side (buckets and unit index updated) then on the row side.
function remove_entry!(M::DualSparseMatrix{T}, r::Int, c::Int) where {T}
    c_rows = M.col_rows[c]
    c_vals = M.col_vals[c]
    idx_r = searchsortedfirst(c_rows, r)
    if idx_r <= length(c_rows) && c_rows[idx_r] == r
        was_unit = is_unit(c_vals[idx_r])
        deleteat!(c_rows, idx_r)
        deleteat!(c_vals, idx_r)
        update_degrees!(M.col_buckets, c, length(c_rows))
        ui = M.unit_index
        unit_index_sync!(ui, c, ui.unit_count[c] - (was_unit ? 1 : 0), length(c_rows))
    end

    r_cols = M.row_cols[r]
    r_vals = M.row_vals[r]
    idx_c = searchsortedfirst(r_cols, c)
    if idx_c <= length(r_cols) && r_cols[idx_c] == c
        deleteat!(r_cols, idx_c)
        deleteat!(r_vals, idx_c)
        update_degrees!(M.row_buckets, r, length(r_cols))
    end
    return nothing
end

####################################
### Structured elimination (SGE) ###
####################################

# (No fill-in is created here, hence no overflow risk.)

function eliminate_structured!(M::DualSparseMatrix{T})::Int where {T}
    n_ones = 0
    progress = true

    deg1_col_bidx = degree_to_bucket_idx(M.col_buckets, 1)
    deg1_row_bidx = degree_to_bucket_idx(M.row_buckets, 1)

    while progress
        progress = false

        for c ∈ copy(M.col_buckets.buckets[deg1_col_bidx])
            @inbounds c_rows = M.col_rows[c]
            length(c_rows) == 1 || continue
            @inbounds r = c_rows[1]
            @inbounds val_pivot = M.col_vals[c][1]
            is_unit(val_pivot) || continue

            empty!(M._sister_rows_buf)
            r_cols_source = M.row_cols[r]
            for i ∈ eachindex(r_cols_source)
                @inbounds push!(M._sister_rows_buf, r_cols_source[i])
            end
            for i ∈ eachindex(M._sister_rows_buf)
                @inbounds sc = M._sister_rows_buf[i]
                sc == c && continue
                remove_entry!(M, r, sc)
            end
            remove_entry!(M, r, c)
            n_ones += 1
            progress = true
        end

        for r ∈ copy(M.row_buckets.buckets[deg1_row_bidx])
            @inbounds r_cols = M.row_cols[r]
            length(r_cols) == 1 || continue
            @inbounds c = r_cols[1]
            @inbounds val_pivot = M.row_vals[r][1]
            is_unit(val_pivot) || continue

            empty!(M._sister_rows_buf)
            c_rows_source = M.col_rows[c]
            for i ∈ eachindex(c_rows_source)
                @inbounds push!(M._sister_rows_buf, c_rows_source[i])
            end
            for i ∈ eachindex(M._sister_rows_buf)
                @inbounds sr = M._sister_rows_buf[i]
                sr == r && continue
                remove_entry!(M, sr, c)
            end
            remove_entry!(M, r, c)
            n_ones += 1
            progress = true
        end
    end

    return n_ones
end

###########################################################
### Checked arithmetic and batched sorted O(n+k) merges ###
###########################################################

# Written directly into the shared scratch buffer (no separate compute/commit phases).

# _sub_mul / _neg_mul / _checked_negate: checked arithmetic helpers, for the signed integers
# of fixed width (those of Base, or wider ones such as Int256 of BitIntegers.jl, whose
# *_with_overflow methods Base.checked_* use), plain arithmetic for BigInt.
@inline _sub_mul(curr::T, factor::T, other::T) where {T<:Signed} =
    Base.checked_sub(curr, Base.checked_mul(factor, other))
@inline _sub_mul(curr::BigInt, factor::BigInt, other::BigInt) = curr - factor * other

@inline _neg_mul(factor::T, other::T) where {T<:Signed} =
    Base.checked_sub(zero(T), Base.checked_mul(factor, other))
@inline _neg_mul(factor::BigInt, other::BigInt) = -(factor * other)

@inline _checked_negate(v::T) where {T<:Signed} = Base.checked_sub(zero(T), v)
@inline _checked_negate(v::BigInt) = -v

"""
    _ensure_capacity!(buf, needed)

Ensure `length(buf) >= needed`, with amortized doubling.
"""
@inline function _ensure_capacity!(buf::Vector{T}, needed::Int) where {T}
    if length(buf) < needed
        resize!(buf, max(needed, 2 * length(buf)))
    end
    return nothing
end

@inline function _merge_write_or_skip!(out_keys::Vector{Int}, out_vals::Vector{T},
                                        k::Int, key::Int, val::T)::Int where {T}
    iszero(val) && return k
    k += 1
    @inbounds out_keys[k] = key
    @inbounds out_vals[k] = val
    return k
end

"""
    merge_row!(M, sr, skip_col, new_cols, new_vals, factor, release)

Sorted merge in O(n_old + n_new), written into the shared scratch buffer (`M._merge_cols_buf` / `M._merge_vals_buf`), then copied into `row_cols[sr]` / `row_vals[sr]`. Checked arithmetic (`_sub_mul` / `_neg_mul`). `release`: the setting `release_empty_vectors` of [`ReductionConfig`](@ref) (a row that becomes empty gets fresh vectors, which frees its capacity).

⚠️ NO guarantee of non-mutation on overflow: `M` MUST be abandoned entirely as soon as an `OverflowError` is raised here, never reused — this is what `reduce_to_core!` guarantees.
"""
function merge_row!(M::DualSparseMatrix{T}, sr::Int, skip_col::Int,
                           new_cols::Vector{Int}, new_vals::Vector{T}, factor::T,
                           release::Bool) where {T}
    old_cols = M.row_cols[sr]
    old_vals = M.row_vals[sr]
    n_old = length(old_cols)
    n_new = length(new_cols)
    old_deg = n_old

    out_cols = M._merge_cols_buf
    out_vals = M._merge_vals_buf
    cap_needed = n_old + n_new
    _ensure_capacity!(out_cols, cap_needed)
    _ensure_capacity!(out_vals, cap_needed)
    k = 0

    i = 1; j = 1
    @inbounds while i <= n_old && j <= n_new
        ok = old_cols[i]
        nk = new_cols[j]
        if ok == skip_col
            i += 1
        elseif ok < nk
            k = _merge_write_or_skip!(out_cols, out_vals, k, ok, old_vals[i])
            i += 1
        elseif ok > nk
            val = _neg_mul(factor, new_vals[j])
            k = _merge_write_or_skip!(out_cols, out_vals, k, nk, val)
            j += 1
        else # ok == nk: the new value entirely replaces the old one
            val = _sub_mul(old_vals[i], factor, new_vals[j])
            k = _merge_write_or_skip!(out_cols, out_vals, k, nk, val)
            i += 1; j += 1
        end
    end
    @inbounds while i <= n_old
        ok = old_cols[i]
        ok != skip_col && (k = _merge_write_or_skip!(out_cols, out_vals, k, ok, old_vals[i]))
        i += 1
    end
    @inbounds while j <= n_new
        val = _neg_mul(factor, new_vals[j])
        k = _merge_write_or_skip!(out_cols, out_vals, k, new_cols[j], val)
        j += 1
    end

    if k == 0 && release
        # Row became empty: replace its vectors by fresh ones. `resize!(v, 0)` would keep the whole
        # capacity (ZT2w4m5 d₄: ~15 GiB of capacity in 1.3 M empty vectors, 44 s of GC).
        # Writing to a slot owned by this row: safe inside a parallel region.
        old_cols = Int[]; old_vals = T[]
        M.row_cols[sr] = old_cols; M.row_vals[sr] = old_vals
    else
        resize!(old_cols, k); copyto!(old_cols, 1, out_cols, 1, k)
        resize!(old_vals, k); copyto!(old_vals, 1, out_vals, 1, k)
    end

    new_deg = length(old_cols)
    new_deg != old_deg && update_degrees!(M.row_buckets, sr, new_deg)
    return nothing
end

"""
    merge_col!(M, c_piv, skip_row, new_rows, factors, v_p, release)

Column-side counterpart of `merge_row!`.
"""
function merge_col!(M::DualSparseMatrix{T}, c_piv::Int, skip_row::Int,
                           new_rows::Vector{Int}, factors::Vector{T}, v_p::T,
                           release::Bool) where {T}
    old_rows = M.col_rows[c_piv]
    old_vals = M.col_vals[c_piv]
    n_old = length(old_rows)
    n_new = length(new_rows)
    old_deg = n_old

    out_rows = M._merge_cols_buf
    out_vals = M._merge_vals_buf
    cap_needed = n_old + n_new
    _ensure_capacity!(out_rows, cap_needed)
    _ensure_capacity!(out_vals, cap_needed)
    k = 0

    i = 1; j = 1
    @inbounds while i <= n_old && j <= n_new
        ok = old_rows[i]
        nk = new_rows[j]
        if ok == skip_row
            i += 1
        elseif ok < nk
            k = _merge_write_or_skip!(out_rows, out_vals, k, ok, old_vals[i])
            i += 1
        elseif ok > nk
            val = _neg_mul(factors[j], v_p)
            k = _merge_write_or_skip!(out_rows, out_vals, k, nk, val)
            j += 1
        else
            val = _sub_mul(old_vals[i], factors[j], v_p)
            k = _merge_write_or_skip!(out_rows, out_vals, k, nk, val)
            i += 1; j += 1
        end
    end
    @inbounds while i <= n_old
        ok = old_rows[i]
        ok != skip_row && (k = _merge_write_or_skip!(out_rows, out_vals, k, ok, old_vals[i]))
        i += 1
    end
    @inbounds while j <= n_new
        val = _neg_mul(factors[j], v_p)
        k = _merge_write_or_skip!(out_rows, out_vals, k, new_rows[j], val)
        j += 1
    end

    if k == 0 && release
        # Column became empty: same replacement as for rows (see `merge_row!`).
        old_rows = Int[]; old_vals = T[]
        M.col_rows[c_piv] = old_rows; M.col_vals[c_piv] = old_vals
    else
        resize!(old_rows, k); copyto!(old_rows, 1, out_rows, 1, k)
        resize!(old_vals, k); copyto!(old_vals, 1, out_vals, 1, k)
    end

    new_deg = length(old_rows)
    new_deg != old_deg && update_degrees!(M.col_buckets, c_piv, new_deg)

    new_unit_count = 0
    @inbounds for v ∈ old_vals
        is_unit(v) && (new_unit_count += 1)
    end
    unit_index_sync!(M.unit_index, c_piv, new_unit_count, new_deg)

    return nothing
end

##########################################################
### Configuration, optional instrumentation, workspace ###
##########################################################

"""
    SparseElementaryDivisors.ReductionTrace()

Trajectory fingerprint, a developer tool (public, not exported): passed to
[`sparse_elementary_divisors`](@ref) with the keyword `trace`, it receives a running hash
(FNV-1a on 64-bit words) of the Markowitz pivots `(row, column)` and of the calls of the
structured elimination, in order. Two computations with the same fingerprint chose the same
pivots in the same order; the test suite compares it with references computed by v0.1.0, for
1 to 4 workers and with the parallel paths forced. Reset at the start of each attempt of the
widening: it describes the attempt retained.

Fields: `h` (hash), `n_pivots` (Markowitz pivots), `n_sge` (calls of the structured
elimination), `sum_ones` (unit pivots found by the structured elimination).
"""
mutable struct ReductionTrace
    h::UInt64
    n_pivots::Int
    n_sge::Int
    sum_ones::Int
end
ReductionTrace() = ReductionTrace(0xcbf29ce484222325, 0, 0, 0)

# Back to the initial state (at the start of each attempt of the widening).
function _reset!(tr::ReductionTrace)
    tr.h = 0xcbf29ce484222325
    tr.n_pivots = 0
    tr.n_sge = 0
    tr.sum_ones = 0
    return tr
end

@inline function _trace_mix!(tr::ReductionTrace, x::Int)
    h = (tr.h ⊻ (x % UInt64)) * 0x00000100000001b3
    tr.h = h ⊻ (h >> 29)
    return nothing
end
@inline function _trace_pivot!(tr::ReductionTrace, p_r::Int, p_c::Int)
    _trace_mix!(tr, p_r); _trace_mix!(tr, p_c)
    tr.n_pivots += 1
    return nothing
end
@inline function _trace_sge!(tr::ReductionTrace, n_ones::Int)
    _trace_mix!(tr, -n_ones - 1)
    tr.n_sge += 1
    tr.sum_ones += n_ones
    return nothing
end
Base.show(io::IO, tr::ReductionTrace) =
    print(io, "ReductionTrace(h=0x", string(tr.h; base = 16), ", pivots=", tr.n_pivots,
          ", sge=", tr.n_sge, ", ones_sge=", tr.sum_ones, ")")

"""
    ReductionStats()

Statistics of a computation, filled in place by [`sparse_elementary_divisors`](@ref) through its
keyword `stats`:

```julia
st = ReductionStats()
d = sparse_elementary_divisors(A; stats = st)
st.value_type        # the integer type retained by the widening
```

All the fields are reset at the start of each call. The counters, times and core describe the
attempt retained (the one that did not overflow); with `value_type = :auto` or a sequence of
types, the attempts that overflowed before it are listed in `failed_attempts`. One object per
computation: two computations running at the same time must not share it.

Widening:
- `value_type`: the type retained, `nothing` before the end of a computation (and after an
  error);
- `failed_attempts`: the attempts that overflowed, as `(value_type, unit_pivots, seconds)`, in
  order; the last one included when the last type does not suffice (the error is then
  raised).

Elimination (attempt retained):
- `unit_pivots`: ±1 pivots of the structured elimination and of the Markowitz loop, that is
  the elementary divisors equal to 1 found before the core;
- `markowitz_pivots`: those of the Markowitz loop;
- `parallel_searches`, `sequential_searches`: bucket scans of the pivot search, parallel or
  sequential;
- `parallel_applications`, `sequential_applications`: pivot applications (merges), parallel
  or sequential.

Times, in seconds:
- `sge_seconds`, `search_seconds`, `apply_seconds`: structured elimination, pivot search,
  pivot application;
- `core_seconds`: extraction of the core and its elementary divisors by Hecke.

Core and context:
- `core_size`: `(rows, cols, nnz)` of the core left by the elimination;
- `core_format`: `:none` (empty core), `:sparse` or `:dense` (see `core_format` in
  [`ReductionConfig`](@ref));
- `nworkers`: number of workers used.
"""
mutable struct ReductionStats
    value_type::Union{Nothing, Type}
    failed_attempts::Vector{NamedTuple{(:value_type, :unit_pivots, :seconds), Tuple{Type, Int, Float64}}}
    unit_pivots::Int
    markowitz_pivots::Int
    parallel_searches::Int
    sequential_searches::Int
    parallel_applications::Int
    sequential_applications::Int
    sge_seconds::Float64
    search_seconds::Float64
    apply_seconds::Float64
    core_seconds::Float64
    core_size::NamedTuple{(:rows, :cols, :nnz), Tuple{Int, Int, Int}}
    core_format::Symbol
    nworkers::Int
end
ReductionStats() =
    ReductionStats(nothing, NamedTuple{(:value_type, :unit_pivots, :seconds), Tuple{Type, Int, Float64}}[],
                   0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, (rows = 0, cols = 0, nnz = 0), :none, 0)

# The counters, times, core and workers back to zero (at the start of each attempt).
function _reset_counters!(s::ReductionStats)
    s.unit_pivots = s.markowitz_pivots = 0
    s.parallel_searches = s.sequential_searches = 0
    s.parallel_applications = s.sequential_applications = 0
    s.sge_seconds = s.search_seconds = s.apply_seconds = s.core_seconds = 0.0
    s.core_size = (rows = 0, cols = 0, nnz = 0)
    s.core_format = :none
    s.nworkers = 0
    return s
end

# Everything back to the state of ReductionStats() (at the start of each call).
function _reset!(s::ReductionStats)
    s.value_type = nothing
    empty!(s.failed_attempts)
    return _reset_counters!(s)
end

_core_string(s::ReductionStats) = string(s.core_size.rows, "×", s.core_size.cols)

function Base.show(io::IO, s::ReductionStats)
    print(io, "ReductionStats(value_type = ", s.value_type === nothing ? "nothing" : s.value_type,
          ", unit_pivots = ", s.unit_pivots, ", core = ", _core_string(s))
    isempty(s.failed_attempts) ||
        print(io, ", failed = [", join((a.value_type for a in s.failed_attempts), ", "), "]")
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", s::ReductionStats)
    r3(x) = round(x; digits = 3)
    println(io, "ReductionStats")
    print(io, "  value type     ", s.value_type === nothing ? "nothing" : s.value_type)
    for a in s.failed_attempts
        print(io, "\n  failed         ", a.value_type, " after ", a.unit_pivots, " unit pivots (",
              r3(a.seconds), " s)")
    end
    print(io, "\n  unit pivots    ", s.unit_pivots, " (Markowitz ", s.markowitz_pivots, ")")
    print(io, "\n  bucket scans   parallel ", s.parallel_searches, ", sequential ", s.sequential_searches)
    print(io, "\n  applications   parallel ", s.parallel_applications, ", sequential ",
          s.sequential_applications)
    print(io, "\n  times (s)      structured ", r3(s.sge_seconds), ", search ", r3(s.search_seconds),
          ", application ", r3(s.apply_seconds), ", core ", r3(s.core_seconds))
    print(io, "\n  core           ", _core_string(s), ", ", s.core_size.nnz, " entries, ", s.core_format)
    print(io, "\n  workers        ", s.nworkers)
end

"""
    ReductionConfig(; kwargs...)

Settings of one attempt of the reduction, each with its default; an immutable value, which can
be kept and reused (for instance for all the differentials of a complex). Pass it to
[`sparse_elementary_divisors`](@ref) with the keyword `config`; each field can also be given
there directly as a keyword, which then takes precedence over `config`:

```julia
cfg = ReductionConfig(nworkers = 4, core_format = :sparse)
sparse_elementary_divisors(A; config = cfg)
sparse_elementary_divisors(A; config = cfg, nworkers = 1)     # nworkers = 1 wins
sparse_elementary_divisors(A; search_mode = :never)           # without building a config
```

No setting changes the elementary divisors; the pivot search parameters change the order of
the pivots (hence the time, the fill-in and the size of the coefficients), the others only how
the work is done. The values are checked at construction (`ArgumentError`).

Pivot search (±1 pivots, Markowitz score `(row degree - 1)(column degree - 1)`, columns visited
by increasing degree):
- `acceptable_score_threshold = 4`: the search stops at the first pivot whose score is at
  most this value (`≥ 0`);
- `max_buckets_scanned = 1`: maximum number of non-empty degree buckets visited by a search
  (`≥ 1`).

Parallelism (the result does not depend on it, bit for bit):
- `nworkers = nothing`: number of workers, the calling thread included; `nothing` means
  `Threads.nthreads()` at the time of the call (start Julia with `-t auto` or `-t N`);
- `apply_min_work = 25_000`: work (entries traversed by the merges of a pivot) below which the
  application of a pivot stays sequential;
- `apply_chunk_max = 32`: maximum size of an application chunk;
- `search_mode = :auto`: `:auto` (cost rule below), `:always` (parallel as soon as
  `nworkers > 1`), `:never`, or `:floor` (parallel when the bucket has at least
  `search_parallel_min_work` entries, without the cost rule);
- `search_parallel_min_work = 5000`: floor of the `:auto` rule, in bucket entries;
- `search_ns_per_entry = 8.0`, `parallel_fixed_cost_ns = 50_000.0`,
  `apply_parallel_penalty = 0.1`: constants of the `:auto` rule, calibrated on N8 (24 threads):
  a bucket of `e` entries is scanned in parallel when the saving
  `e · search_ns_per_entry · (1 - 1/nworkers)` exceeds
  `parallel_fixed_cost_ns + apply_parallel_penalty · (recent application time in ns)`;
- `apply_ema_alpha = 0.05`: weight of the moving average of the application time, in
  `(0, 1]`;
- `search_chunk_min = 4`, `search_chunk_max = 256`: bounds on the size of the search chunks
  (columns).

Core (what the elimination leaves, passed to Hecke):
- `core_format = :auto`: `:sparse` (`SMat`), `:dense` (Flint `ZZMatrix`, an extraction without
  pressure on the garbage collector: ZT2w4m5 d₄, 0.6 s instead of 422 s), or `:auto` (dense
  when the density of the core is at least `dense_min_density` and its dense form takes at most
  `dense_max_bytes` bytes, sparse otherwise);
- `dense_min_density = 0.05`, in `[0, 1]`;
- `dense_max_bytes = 16 * 2^30`.

Memory:
- `release_empty_vectors = false`: if `true`, a row or column that becomes empty during a merge
  gets fresh vectors, which frees its capacity (by default its vectors are emptied and keep their
  capacity). On Windows, this reduced the heap of ZT2w4m5 d₄ from 22 to 9 GiB but made the end
  of the reduction twice as slow (garbage collection 443 to 460 s instead of about 44 s; total
  about 860 to 890 s instead of 434 s); not measured on Linux yet.
"""
struct ReductionConfig
    acceptable_score_threshold::Int
    max_buckets_scanned::Int
    nworkers::Union{Nothing, Int}
    apply_min_work::Int
    apply_chunk_max::Int
    search_mode::Symbol
    search_parallel_min_work::Int
    search_ns_per_entry::Float64
    parallel_fixed_cost_ns::Float64
    apply_parallel_penalty::Float64
    apply_ema_alpha::Float64
    search_chunk_min::Int
    search_chunk_max::Int
    core_format::Symbol
    dense_min_density::Float64
    dense_max_bytes::Int
    release_empty_vectors::Bool

    # The only constructor: every other one ends here, so that every configuration is checked.
    function ReductionConfig(acceptable_score_threshold, max_buckets_scanned, nworkers,
                             apply_min_work, apply_chunk_max, search_mode, search_parallel_min_work,
                             search_ns_per_entry, parallel_fixed_cost_ns, apply_parallel_penalty,
                             apply_ema_alpha, search_chunk_min, search_chunk_max, core_format,
                             dense_min_density, dense_max_bytes, release_empty_vectors)
        c = new(acceptable_score_threshold, max_buckets_scanned, nworkers, apply_min_work,
                apply_chunk_max, search_mode, search_parallel_min_work, search_ns_per_entry,
                parallel_fixed_cost_ns, apply_parallel_penalty, apply_ema_alpha, search_chunk_min,
                search_chunk_max, core_format, dense_min_density, dense_max_bytes,
                release_empty_vectors)
        _check_config(c)
        return c
    end
end

function ReductionConfig(; acceptable_score_threshold = 4, max_buckets_scanned = 1,
                         nworkers = nothing, apply_min_work = 25_000, apply_chunk_max = 32,
                         search_mode = :auto, search_parallel_min_work = 5000,
                         search_ns_per_entry = 8.0, parallel_fixed_cost_ns = 50_000.0,
                         apply_parallel_penalty = 0.1, apply_ema_alpha = 0.05,
                         search_chunk_min = 4, search_chunk_max = 256, core_format = :auto,
                         dense_min_density = 0.05, dense_max_bytes = 16 * 2^30,
                         release_empty_vectors = false)
    return ReductionConfig(acceptable_score_threshold, max_buckets_scanned, nworkers,
                           apply_min_work, apply_chunk_max, search_mode, search_parallel_min_work,
                           search_ns_per_entry, parallel_fixed_cost_ns, apply_parallel_penalty,
                           apply_ema_alpha, search_chunk_min, search_chunk_max, core_format,
                           dense_min_density, dense_max_bytes, release_empty_vectors)
end

function _check_config(c::ReductionConfig)
    bad(field, rule) = throw(ArgumentError("ReductionConfig: $field must be $rule (got " *
                                           "$(repr(getfield(c, field))))"))
    c.acceptable_score_threshold >= 0 || bad(:acceptable_score_threshold, "≥ 0")
    c.max_buckets_scanned >= 1 || bad(:max_buckets_scanned, "≥ 1")
    (c.nworkers === nothing || c.nworkers >= 1) || bad(:nworkers, "nothing or ≥ 1")
    c.apply_min_work >= 0 || bad(:apply_min_work, "≥ 0")
    c.apply_chunk_max >= 1 || bad(:apply_chunk_max, "≥ 1")
    c.search_mode ∈ (:auto, :always, :never, :floor) ||
        bad(:search_mode, ":auto, :always, :never or :floor")
    c.search_parallel_min_work >= 0 || bad(:search_parallel_min_work, "≥ 0")
    c.search_ns_per_entry >= 0 || bad(:search_ns_per_entry, "≥ 0")
    c.parallel_fixed_cost_ns >= 0 || bad(:parallel_fixed_cost_ns, "≥ 0")
    c.apply_parallel_penalty >= 0 || bad(:apply_parallel_penalty, "≥ 0")
    0 < c.apply_ema_alpha <= 1 || bad(:apply_ema_alpha, "in (0, 1]")
    c.search_chunk_min >= 1 || bad(:search_chunk_min, "≥ 1")
    c.search_chunk_max >= c.search_chunk_min || bad(:search_chunk_max, "≥ search_chunk_min")
    c.core_format ∈ (:auto, :sparse, :dense) || bad(:core_format, ":auto, :sparse or :dense")
    0 <= c.dense_min_density <= 1 || bad(:dense_min_density, "in [0, 1]")
    c.dense_max_bytes >= 0 || bad(:dense_max_bytes, "≥ 0")
    return nothing
end

# Only the fields that differ from the default, in the form of the keyword constructor.
function Base.show(io::IO, c::ReductionConfig)
    d = ReductionConfig()
    changed = [f for f in fieldnames(ReductionConfig) if getfield(c, f) != getfield(d, f)]
    print(io, "ReductionConfig(")
    join(io, (string(f, " = ", repr(getfield(c, f))) for f in changed), ", ")
    print(io, ")")
end

# Number of workers of a configuration, `nothing` resolved at the time of the call.
_nworkers(c::ReductionConfig) = c.nworkers === nothing ? Threads.nthreads() : c.nworkers

# `config` with the keywords of `settings` (fields of ReductionConfig) replacing its fields.
function _with_settings(config::ReductionConfig, settings)
    isempty(settings) && return config
    names = fieldnames(ReductionConfig)
    for k in keys(settings)
        k ∈ names || throw(ArgumentError(
            "unknown keyword argument $k; the keywords are value_type, config, stats, trace, " *
            "probe, verbose and the fields of ReductionConfig: $(join(names, ", "))"))
    end
    vals = map(f -> haskey(settings, f) ? settings[f] : getfield(config, f), names)
    return ReductionConfig(vals...)
end

"Per-chunk results of a parallel search."
mutable struct SearchWorkspace
    cap::Int
    status::Vector{Int8}      # 0 = chunk not processed, 1 = processed
    good::Vector{Bool}        # the chunk reached the threshold (min_* = its FIRST candidate <= threshold)
    min_score::Vector{Int}
    min_r::Vector{Int}
    min_c::Vector{Int}
end
SearchWorkspace(cap::Int = 16) = SearchWorkspace(cap, zeros(Int8, cap), fill(false, cap),
                                   fill(typemax(Int), cap), zeros(Int, cap), zeros(Int, cap))

function _sws_ensure!(sws::SearchWorkspace, n::Int)
    n <= sws.cap && return nothing
    newcap = max(n, 2 * sws.cap)
    resize!(sws.status, newcap); resize!(sws.good, newcap)
    resize!(sws.min_score, newcap); resize!(sws.min_r, newcap); resize!(sws.min_c, newcap)
    sws.cap = newcap
    return nothing
end

"""
Workspace of a reduction: per-worker merge buffers, unit-count buffer, search results, EMA, and
the times accumulated for `ReductionStats` (integer nanoseconds, converted to seconds at the end
of `reduce_to_core!`).
"""
mutable struct ReductionWorkspace{T}
    nw::Int
    wkeys::Vector{Vector{Int}}
    wvals::Vector{Vector{T}}
    col_unit_buf::Vector{Int}
    sws::SearchWorkspace
    apply_ema::Float64        # recent application duration (ns), exponential moving average
    ema_init::Bool
    sge_ns::Int
    search_ns::Int
    apply_ns::Int
end
function ReductionWorkspace{T}(cfg::ReductionConfig) where {T}
    nw = _nworkers(cfg)
    return ReductionWorkspace{T}(nw, [Int[] for _ ∈ 1:nw], [T[] for _ ∈ 1:nw], Int[], SearchWorkspace(),
                                 0.0, false, 0, 0, 0)
end

@inline function _update_ema!(wsp::ReductionWorkspace, dt_ns::Int, alpha::Float64)
    if wsp.ema_init
        wsp.apply_ema += alpha * (dt_ns - wsp.apply_ema)
    else
        wsp.apply_ema = Float64(dt_ns); wsp.ema_init = true
    end
    return nothing
end

######################################################
### Self-scheduling parallel loop (atomic counter) ###
######################################################

"""
    _parallel_dynamic!(body, n, nw, chunk)

Run `body(t, w)` for t ∈ 1:n on `nw` workers (w ∈ 1:nw; the calling thread is worker 1), dynamic distribution by chunks of `chunk` elements (atomic counter). `w` is assigned HERE (stable for the task): safe to index private buffers even if tasks migrate.

Exception robustness: each worker catches the exception, stores it (only the first) and raises a stop flag checked by the others between two chunks; ALL tasks are awaited, then the original exception is rethrown by the caller (never a `TaskFailedException`).
"""
function _parallel_dynamic!(body::F, n::Int, nw::Int, chunk::Int) where {F}
    counter = Threads.Atomic{Int}(0)
    abort = Threads.Atomic{Bool}(false)
    err = Ref{Any}(nothing)
    err_lock = ReentrantLock()
    worker = function (w::Int)
        try
            while !abort[]
                lo = Threads.atomic_add!(counter, chunk)
                lo >= n && break
                hi = min(lo + chunk, n)
                for t ∈ (lo + 1):hi
                    body(t, w)
                end
            end
        catch e
            lock(err_lock) do
                err[] === nothing && (err[] = e)
            end
            abort[] = true
        end
        return nothing
    end
    tasks = Vector{Task}(undef, nw - 1)
    for w ∈ 2:nw
        let ww = w
            tasks[w - 1] = Threads.@spawn worker(ww)
        end
    end
    worker(1)
    for tk ∈ tasks
        wait(tk)
    end
    err[] === nothing || throw(err[])
    return nothing
end

#####################################################################
### Per-worker merges (explicit output buffers, NO bucket update) ###
#####################################################################

function _merge_row_buffered!(M::DualSparseMatrix{T}, sr::Int, skip_col::Int,
                       new_cols::Vector{Int}, new_vals::Vector{T}, factor::T,
                       out_cols::Vector{Int}, out_vals::Vector{T}, release::Bool) where {T}
    old_cols = M.row_cols[sr]
    old_vals = M.row_vals[sr]
    n_old = length(old_cols)
    n_new = length(new_cols)
    cap_needed = n_old + n_new
    _ensure_capacity!(out_cols, cap_needed)
    _ensure_capacity!(out_vals, cap_needed)
    k = 0
    i = 1; j = 1
    @inbounds while i <= n_old && j <= n_new
        ok = old_cols[i]
        nk = new_cols[j]
        if ok == skip_col
            i += 1
        elseif ok < nk
            k = _merge_write_or_skip!(out_cols, out_vals, k, ok, old_vals[i])
            i += 1
        elseif ok > nk
            val = _neg_mul(factor, new_vals[j])
            k = _merge_write_or_skip!(out_cols, out_vals, k, nk, val)
            j += 1
        else
            val = _sub_mul(old_vals[i], factor, new_vals[j])
            k = _merge_write_or_skip!(out_cols, out_vals, k, nk, val)
            i += 1; j += 1
        end
    end
    @inbounds while i <= n_old
        ok = old_cols[i]
        ok != skip_col && (k = _merge_write_or_skip!(out_cols, out_vals, k, ok, old_vals[i]))
        i += 1
    end
    @inbounds while j <= n_new
        val = _neg_mul(factor, new_vals[j])
        k = _merge_write_or_skip!(out_cols, out_vals, k, new_cols[j], val)
        j += 1
    end
    if k == 0 && release
        # Row became empty: replace its vectors by fresh ones. `resize!(v, 0)` would keep the whole
        # capacity (ZT2w4m5 d₄: ~15 GiB of capacity in 1.3 M empty vectors, 44 s of GC).
        # Writing to a slot owned by this row: safe inside a parallel region.
        old_cols = Int[]; old_vals = T[]
        M.row_cols[sr] = old_cols; M.row_vals[sr] = old_vals
    else
        resize!(old_cols, k); copyto!(old_cols, 1, out_cols, 1, k)
        resize!(old_vals, k); copyto!(old_vals, 1, out_vals, 1, k)
    end
    return nothing
end

# Returns the number of ±1 entries of the column after the merge.
function _merge_col_buffered!(M::DualSparseMatrix{T}, c_piv::Int, skip_row::Int,
                       new_rows::Vector{Int}, factors::Vector{T}, v_p::T,
                       out_rows::Vector{Int}, out_vals::Vector{T}, release::Bool)::Int where {T}
    old_rows = M.col_rows[c_piv]
    old_vals = M.col_vals[c_piv]
    n_old = length(old_rows)
    n_new = length(new_rows)
    cap_needed = n_old + n_new
    _ensure_capacity!(out_rows, cap_needed)
    _ensure_capacity!(out_vals, cap_needed)
    k = 0
    i = 1; j = 1
    @inbounds while i <= n_old && j <= n_new
        ok = old_rows[i]
        nk = new_rows[j]
        if ok == skip_row
            i += 1
        elseif ok < nk
            k = _merge_write_or_skip!(out_rows, out_vals, k, ok, old_vals[i])
            i += 1
        elseif ok > nk
            val = _neg_mul(factors[j], v_p)
            k = _merge_write_or_skip!(out_rows, out_vals, k, nk, val)
            j += 1
        else
            val = _sub_mul(old_vals[i], factors[j], v_p)
            k = _merge_write_or_skip!(out_rows, out_vals, k, nk, val)
            i += 1; j += 1
        end
    end
    @inbounds while i <= n_old
        ok = old_rows[i]
        ok != skip_row && (k = _merge_write_or_skip!(out_rows, out_vals, k, ok, old_vals[i]))
        i += 1
    end
    @inbounds while j <= n_new
        val = _neg_mul(factors[j], v_p)
        k = _merge_write_or_skip!(out_rows, out_vals, k, new_rows[j], val)
        j += 1
    end
    if k == 0 && release
        # Column became empty: same replacement as for rows (see `merge_row!`).
        old_rows = Int[]; old_vals = T[]
        M.col_rows[c_piv] = old_rows; M.col_vals[c_piv] = old_vals
    else
        resize!(old_rows, k); copyto!(old_rows, 1, out_rows, 1, k)
        resize!(old_vals, k); copyto!(old_vals, 1, out_vals, 1, k)
    end
    new_unit_count = 0
    @inbounds for v ∈ old_vals
        is_unit(v) && (new_unit_count += 1)
    end
    return new_unit_count
end

########################################################################
### Markowitz pivoting: search (sequential / parallel) + application ###
########################################################################

# --------------------------- #
# Parallel search of a bucket #
# --------------------------- #

# Sequential scan of chunk `ch` (positions k_lo:k_hi of the bucket).
function _scan_chunk!(M::DualSparseMatrix{T}, bucket::Vector{Int}, target_deg::Int, thr::Int,
                      ch::Int, chunk::Int, sws::SearchWorkspace) where {T}
    nb::Int = length(bucket)
    k_lo::Int = (ch - 1) * chunk + 1
    k_hi::Int = min(ch * chunk, nb)
    best_s::Int = typemax(Int)
    best_r::Int = -1
    best_c::Int = -1
    good::Bool = false
    for k::Int ∈ k_lo:k_hi
        @inbounds c::Int = bucket[k]
        @inbounds deg_c::Int = length(M.col_rows[c])
        deg_c != target_deg && continue
        c_rows::Vector{Int} = M.col_rows[c]
        c_vals::Vector{T} = M.col_vals[c]
        for idx_r::Int ∈ 1:deg_c
            @inbounds val::T = c_vals[idx_r]
            if is_unit(val)
                @inbounds r::Int = c_rows[idx_r]
                @inbounds deg_r::Int = length(M.row_cols[r])
                score::Int = (deg_r - 1) * (deg_c - 1)
                if score < best_s
                    best_s = score
                    best_r = r
                    best_c = c
                end
                if best_s <= thr
                    good = true
                    break
                end
            end
        end
        good && break
    end
    @inbounds begin
        sws.good[ch] = good
        sws.min_score[ch] = best_s
        sws.min_r[ch] = best_r
        sws.min_c[ch] = best_c
        sws.status[ch] = Int8(1)
    end
    return nothing
end

# Parallel region: returns the number of chunks (results in `sws`).
function _par_scan_bucket!(M::DualSparseMatrix{T}, bucket::Vector{Int}, target_deg::Int, thr::Int,
                           nw::Int, sws::SearchWorkspace, chunk_min::Int, chunk_max::Int)::Int where {T}
    nb::Int = length(bucket)
    chunk::Int = clamp(nb ÷ (8 * nw), chunk_min, chunk_max)
    nch::Int = cld(nb, chunk)
    _sws_ensure!(sws, nch)
    @inbounds for i ∈ 1:nch
        sws.status[i] = Int8(0)
    end
    good_min = Threads.Atomic{Int}(typemax(Int))
    body = function (ch::Int, w::Int)
        ch > good_min[] && return nothing      # an earlier chunk has already reached the threshold
        _scan_chunk!(M, bucket, target_deg, thr, ch, chunk, sws)
        @inbounds if sws.good[ch]
            Threads.atomic_min!(good_min, ch)
        end
        return nothing
    end
    _parallel_dynamic!(body, nch, nw, 1)
    return nch
end

@inline function _use_parallel_search(cfg::ReductionConfig, wsp::ReductionWorkspace, nw::Int, entries::Int)::Bool
    nw <= 1 && return false
    mode = cfg.search_mode
    mode === :never && return false
    mode === :always && return true
    mode === :floor && return entries >= cfg.search_parallel_min_work
    entries < cfg.search_parallel_min_work && return false
    gain = entries * cfg.search_ns_per_entry * (1.0 - 1.0 / nw)
    return gain > cfg.parallel_fixed_cost_ns + cfg.apply_parallel_penalty * wsp.apply_ema
end

# Full search (main buckets, then overflow bucket). Returns (p_r, p_c, best_score).
# Only the scan of a single bucket may be parallel.
function _search_pivot!(M::DualSparseMatrix{T}, wsp::ReductionWorkspace{T}, cfg::ReductionConfig, nw::Int,
                        stats::Union{Nothing, ReductionStats}) where {T}
    acceptable_score_threshold::Int = cfg.acceptable_score_threshold
    max_buckets_scanned::Int = cfg.max_buckets_scanned
    best_score::Int = typemax(Int)
    p_r::Int = -1
    p_c::Int = -1
    found_good_pivot::Bool = false
    n_buckets_scanned::Int = 0

    max_deg::Int = M.col_buckets.max_deg
    unit_buckets::DegreeBuckets = M.unit_index.buckets

    # Cursor (P5): the buckets of the degrees 2 to low - 1 are empty; skip the empty ones that
    # follow, and keep the position for the next search. The non-empty buckets visited are
    # those of the full loop, in the same order: same pivot.
    first_deg::Int = unit_buckets.low
    @inbounds while first_deg <= max_deg && isempty(unit_buckets.buckets[first_deg + 1])
        first_deg += 1
    end
    unit_buckets.low = first_deg

    for target_deg::Int ∈ first_deg:max_deg
        bidx::Int = degree_to_bucket_idx(unit_buckets, target_deg)
        bucket = unit_buckets.buckets[bidx]
        isempty(bucket) && continue
        n_buckets_scanned += 1
        if _use_parallel_search(cfg, wsp, nw, length(bucket) * target_deg)
            stats === nothing || (stats.parallel_searches += 1)
            sws = wsp.sws
            nch::Int = _par_scan_bucket!(M, bucket, target_deg, acceptable_score_threshold, nw, sws,
                                         cfg.search_chunk_min, cfg.search_chunk_max)
            for ch::Int ∈ 1:nch
                @inbounds st = sws.status[ch]
                st == Int8(0) && break
                @inbounds gd = sws.good[ch]
                @inbounds ms = sws.min_score[ch]
                if gd || ms < best_score
                    best_score = ms
                    @inbounds p_r = sws.min_r[ch]
                    @inbounds p_c = sws.min_c[ch]
                end
                if gd
                    found_good_pivot = true
                    break
                end
            end
        else
            stats === nothing || (stats.sequential_searches += 1)
            for c::Int ∈ bucket
                @inbounds deg_c::Int = length(M.col_rows[c])
                deg_c != target_deg && continue
                c_rows::Vector{Int} = M.col_rows[c]
                c_vals::Vector{T} = M.col_vals[c]
                len_c_rows::Int = length(c_rows)
                for idx_r::Int ∈ 1:len_c_rows
                    @inbounds val::T = c_vals[idx_r]
                    if is_unit(val)
                        @inbounds r::Int = c_rows[idx_r]
                        @inbounds deg_r::Int = length(M.row_cols[r])
                        score::Int = (deg_r - 1) * (deg_c - 1)
                        if score < best_score
                            best_score = score
                            p_r = r
                            p_c = c
                        end
                        if best_score <= acceptable_score_threshold
                            found_good_pivot = true
                            break
                        end
                    end
                end
                found_good_pivot && break
            end
        end
        found_good_pivot && break
        n_buckets_scanned >= max_buckets_scanned && break
    end

    # Overflow bucket ("degree > max_deg"): sequential.
    if p_r == -1
        overflow_bidx::Int = degree_to_bucket_idx(unit_buckets, max_deg + 1)
        found_good_pivot_sec::Bool = false
        for c_secours::Int ∈ unit_buckets.buckets[overflow_bidx]
            @inbounds deg_c_sec::Int = length(M.col_rows[c_secours])
            deg_c_sec <= 1 && continue
            c_rows_sec::Vector{Int} = M.col_rows[c_secours]
            c_vals_sec::Vector{T} = M.col_vals[c_secours]
            len_c_rows_sec::Int = length(c_rows_sec)
            for idx_r_sec::Int ∈ 1:len_c_rows_sec
                @inbounds val_sec::T = c_vals_sec[idx_r_sec]
                if is_unit(val_sec)
                    @inbounds r_sec::Int = c_rows_sec[idx_r_sec]
                    @inbounds deg_r_sec::Int = length(M.row_cols[r_sec])
                    score_sec::Int = (deg_r_sec - 1) * (deg_c_sec - 1)
                    if score_sec < best_score
                        best_score = score_sec
                        p_r = r_sec
                        p_c = c_secours
                    end
                    if best_score <= acceptable_score_threshold
                        found_good_pivot_sec = true
                        break
                    end
                end
            end
            found_good_pivot_sec && break
        end
    end
    return p_r, p_c, best_score
end

# ------------------- #
# Fill-in application #
# ------------------- #

# Sequential path: merges with a shared buffer, buckets updated on the fly.
function _apply_pivot_sequential!(M::DualSparseMatrix{T}, p_r::Int, p_c::Int, n_sis::Int, n_pc::Int,
                                  release::Bool) where {T}
    @inbounds for i ∈ 1:n_sis
        sr = M._sister_rows_buf[i]
        factor = M._factors_buf[i]
        merge_row!(M, sr, p_c, M._p_row_cols_buf, M._p_row_vals_buf, factor, release)
    end
    @inbounds for j ∈ 1:n_pc
        c_piv = M._p_row_cols_buf[j]
        v_p = M._p_row_vals_buf[j]
        merge_col!(M, c_piv, p_r, M._sister_rows_buf, M._factors_buf, v_p, release)
    end
    return nothing
end

@inline function _fillin_item!(M::DualSparseMatrix{T}, wsp::ReductionWorkspace{T}, t::Int, w::Int,
                                p_r::Int, p_c::Int, n_sis::Int, release::Bool) where {T}
    out_keys = wsp.wkeys[w]
    out_vals = wsp.wvals[w]
    if t <= n_sis
        @inbounds sr = M._sister_rows_buf[t]
        @inbounds factor = M._factors_buf[t]
        _merge_row_buffered!(M, sr, p_c, M._p_row_cols_buf, M._p_row_vals_buf, factor, out_keys, out_vals,
                             release)
    else
        j = t - n_sis
        @inbounds c_piv = M._p_row_cols_buf[j]
        @inbounds v_p = M._p_row_vals_buf[j]
        cnt = _merge_col_buffered!(M, c_piv, p_r, M._sister_rows_buf, M._factors_buf, v_p, out_keys, out_vals,
                                   release)
        @inbounds wsp.col_unit_buf[j] = cnt
    end
    return nothing
end

# Parallel path: one region (rows AND columns), then a sequential pass over the buckets in the
# same order as the sequential path (the only place where buckets are modified).
function _apply_pivot_parallel!(M::DualSparseMatrix{T}, wsp::ReductionWorkspace{T}, p_r::Int, p_c::Int,
                             n_sis::Int, n_pc::Int, nw::Int, chunk_max::Int, release::Bool) where {T}
    n_tot = n_sis + n_pc
    length(wsp.col_unit_buf) < n_pc && resize!(wsp.col_unit_buf, n_pc)
    chunk = clamp(n_tot ÷ (8 * nw), 1, chunk_max)
    body = (t::Int, w::Int) -> _fillin_item!(M, wsp, t, w, p_r, p_c, n_sis, release)
    _parallel_dynamic!(body, n_tot, nw, chunk)

    @inbounds for i ∈ 1:n_sis
        sr = M._sister_rows_buf[i]
        update_degrees!(M.row_buckets, sr, length(M.row_cols[sr]))
    end
    @inbounds for j ∈ 1:n_pc
        c_piv = M._p_row_cols_buf[j]
        nd = length(M.col_rows[c_piv])
        update_degrees!(M.col_buckets, c_piv, nd)
        unit_index_sync!(M.unit_index, c_piv, wsp.col_unit_buf[j], nd)
    end
    return nothing
end

"""
    markowitz_pivot!(M; config = ReductionConfig(), workspace = nothing, trace = nothing,
                     stats = nothing)::Bool

Perform one pivot step, with the settings of `config` (search parameters, parallelism,
`release_empty_vectors`). Search and application may be parallel; the result is bit-identical
for any number of workers. `workspace`: to provide in order to avoid recreating it at each call
(done by `reduce_to_core!`); with `stats`, the times are accumulated in it (and converted into
`stats` by `reduce_to_core!`), the counters directly in `stats`. `trace`: trajectory
fingerprint.
"""
function markowitz_pivot!(M::DualSparseMatrix{T};
                          config::ReductionConfig = ReductionConfig(),
                          workspace::Union{Nothing, ReductionWorkspace{T}} = nothing,
                          trace::Union{Nothing, ReductionTrace} = nothing,
                          stats::Union{Nothing, ReductionStats} = nothing)::Bool where {T}
    wsp::ReductionWorkspace{T} = workspace === nothing ? ReductionWorkspace{T}(config) : workspace
    nw::Int = wsp.nw
    timing::Bool = stats !== nothing || (nw > 1 && config.search_mode === :auto)
    release::Bool = config.release_empty_vectors

    t0::UInt64 = stats === nothing ? UInt64(0) : time_ns()
    p_r, p_c, best_score = _search_pivot!(M, wsp, config, nw, stats)
    stats === nothing || (wsp.search_ns += Int(time_ns() - t0))
    if p_r == -1; return false; end
    trace === nothing || _trace_pivot!(trace, p_r, p_c)

    # --- Phase 0: capture the sister rows and the pivot value, and estimate the work.
    # work = Σ_sisters (row length) + Σ_columns (column length) + 2·n_sis·n_pc.
    empty!(M._sister_rows_buf); empty!(M._sister_vals_buf)
    src_sister_rows::Vector{Int} = M.col_rows[p_c]
    src_sister_vals::Vector{T} = M.col_vals[p_c]
    val_pivot::T = zero(T)
    work::Int = 0
    for i::Int ∈ eachindex(src_sister_rows)
        @inbounds rr = src_sister_rows[i]
        @inbounds vv = src_sister_vals[i]
        if rr == p_r
            val_pivot = vv
        else
            push!(M._sister_rows_buf, rr)
            push!(M._sister_vals_buf, vv)
            work += length(M.row_cols[rr])
        end
    end
    n_sis::Int = length(M._sister_rows_buf)

    empty!(M._p_row_cols_buf); empty!(M._p_row_vals_buf)
    src_prow_cols::Vector{Int} = M.row_cols[p_r]
    src_prow_vals::Vector{T} = M.row_vals[p_r]
    for j::Int ∈ eachindex(src_prow_cols)
        @inbounds cc = src_prow_cols[j]
        cc == p_c && continue
        @inbounds push!(M._p_row_cols_buf, cc)
        @inbounds push!(M._p_row_vals_buf, src_prow_vals[j])
        work += length(M.col_rows[cc])
    end
    n_pc::Int = length(M._p_row_cols_buf)
    work += 2 * n_sis * n_pc

    # --- Phase 1: factor for each sister row (±sister_val, never a division).
    resize!(M._factors_buf, n_sis)
    @inbounds for i ∈ 1:n_sis
        M._factors_buf[i] = val_pivot == one(T) ? M._sister_vals_buf[i] :
                             _checked_negate(M._sister_vals_buf[i])
    end

    # --- Phases 2/3: application.
    t1::UInt64 = timing ? time_ns() : UInt64(0)
    if nw > 1 && work >= config.apply_min_work
        _apply_pivot_parallel!(M, wsp, p_r, p_c, n_sis, n_pc, nw, config.apply_chunk_max, release)
        stats === nothing || (stats.parallel_applications += 1)
    else
        _apply_pivot_sequential!(M, p_r, p_c, n_sis, n_pc, release)
        stats === nothing || (stats.sequential_applications += 1)
    end
    if timing
        dt = Int(time_ns() - t1)
        nw > 1 && _update_ema!(wsp, dt, config.apply_ema_alpha)
        stats === nothing || (wsp.apply_ns += dt; stats.markowitz_pivots += 1)
    end

    # --- Cleanup.
    # Pivot row and column: REPLACE them by fresh empty vectors instead of `empty!`, which would
    # keep all the allocated capacity (ZT2w4m5 d₄: 21 GiB of capacity for 2.9 GiB of content,
    # 43 s of GC). No alias exists after Phase 0 (the pivot buffers are copies).
    M.row_cols[p_r] = Int[]; M.row_vals[p_r] = T[]
    M.col_rows[p_c] = Int[]; M.col_vals[p_c] = T[]
    update_degrees!(M.row_buckets, p_r, 0)
    update_degrees!(M.col_buckets, p_c, 0)
    unit_index_clear_column!(M.unit_index, p_c)

    return true
end

#######################################################################
### Core extraction (Hecke SMat), reduction and top-level functions ###
#######################################################################

function extract_core_smat(M::DualSparseMatrix)
    active_rows = [r for r ∈ 1:M.num_rows if !isempty(M.row_cols[r])]
    active_cols = [c for c ∈ 1:M.num_cols if !isempty(M.col_rows[c])]

    if isempty(active_rows) || isempty(active_cols)
        return Hecke.sparse_matrix(Hecke.ZZ)
    end

    col_inverse = zeros(Int, M.num_cols)
    for (new_idx, old_idx) ∈ enumerate(active_cols)
        col_inverse[old_idx] = new_idx
    end

    S = Hecke.sparse_matrix(Hecke.ZZ)
    for old_r ∈ active_rows
        hecke_row = Tuple{Int, Hecke.ZZRingElem}[]
        for idx ∈ eachindex(M.row_cols[old_r])
            @inbounds old_c = M.row_cols[old_r][idx]
            @inbounds val = M.row_vals[old_r][idx]
            new_c = col_inverse[old_c]
            new_c > 0 && push!(hecke_row, (new_c, Hecke.ZZ(_to_zz(val))))
        end
        push!(S, Hecke.sparse_row(Hecke.ZZ, hecke_row))
    end
    return S
end

"Recursively unwrap `TaskFailedException`s (defensive: `_parallel_dynamic!` already rethrows the original exception)."
function _unwrap_task_failure(e)
    while e isa TaskFailedException
        e = e.task.result
    end
    return e
end

"""
    reduce_to_core!(M; config = ReductionConfig(), trace = nothing, stats = nothing,
                    probe = nothing)::Int

Structured elimination + ±1 Markowitz loop, with the settings of `config`; returns the number
of unit pivots. An `OverflowError` (including one coming from a parallel task) becomes a
[`CoefficientOverflowError`](@ref) naming `T`; `M` is then inconsistent and MUST NOT be used any
more (the widening restarts from the input). `stats` receives the counters, the times (in
seconds), `unit_pivots` and the number of workers of this reduction; `trace`, the trajectory
fingerprint; `probe`, a function of no argument, is called after each Markowitz pivot
(diagnostics).
"""
function reduce_to_core!(M::DualSparseMatrix{T};
                         config::ReductionConfig = ReductionConfig(),
                         trace::Union{Nothing, ReductionTrace} = nothing,
                         stats::Union{Nothing, ReductionStats} = nothing,
                         probe = nothing)::Int where {T}
    wsp = ReductionWorkspace{T}(config)
    stats === nothing || (stats.nworkers = wsp.nw)
    total_ones = 0
    while true
        t_sge = stats === nothing ? UInt64(0) : time_ns()
        n_sge = eliminate_structured!(M)
        stats === nothing || (wsp.sge_ns += Int(time_ns() - t_sge))
        total_ones += n_sge
        trace === nothing || _trace_sge!(trace, n_sge)
        local pivot_ok::Bool
        try
            pivot_ok = markowitz_pivot!(M; config, workspace = wsp, trace, stats)
        catch e
            ue = _unwrap_task_failure(e)
            ue isa OverflowError || rethrow()
            throw(CoefficientOverflowError(T, total_ones, nothing))
        end
        pivot_ok || break
        total_ones += 1
        probe === nothing || probe()
    end
    if stats !== nothing
        stats.unit_pivots = total_ones
        stats.sge_seconds = wsp.sge_ns / 1e9
        stats.search_seconds = wsp.search_ns / 1e9
        stats.apply_seconds = wsp.apply_ns / 1e9
    end
    return total_ones
end

"""
    _core_dims(M) -> (active rows, active columns, non-zero entries)
"""
function _core_dims(M::DualSparseMatrix)
    nr = 0; nnz = 0
    for r ∈ 1:M.num_rows
        l = length(M.row_cols[r])
        l > 0 && (nr += 1; nnz += l)
    end
    nc = 0
    for c ∈ 1:M.num_cols
        isempty(M.col_rows[c]) || (nc += 1)
    end
    return nr, nc, nnz
end

"""
    _reduce_and_divisors!(M; config = ReductionConfig(), trace = nothing, stats = nothing,
                          probe = nothing, verbose = false)

Reduction (structured elimination + ±1 pivots), then elementary divisors of the core by Hecke,
in the format given by `config.core_format`: `:sparse` (`SMat`), `:dense` (Flint `ZZMatrix`),
`:auto` (dense if the density of the core is at least `config.dense_min_density` and
`8 × rows × columns ≤ config.dense_max_bytes`, sparse otherwise). The result is the same in all
cases. `stats` also receives the size, format and time of the core; `verbose = true` prints the
size and format of the core.
"""
function _reduce_and_divisors!(M::DualSparseMatrix;
                               config::ReductionConfig = ReductionConfig(),
                               trace::Union{Nothing, ReductionTrace} = nothing,
                               stats::Union{Nothing, ReductionStats} = nothing,
                               probe = nothing, verbose::Bool = false)
    n_ones = reduce_to_core!(M; config, trace, stats, probe)

    nr, nc, nnz = _core_dims(M)
    stats === nothing || (stats.core_size = (rows = nr, cols = nc, nnz = nnz))
    verbose && println("Core size: $(nr)×$(nc)")
    if nr == 0 || nc == 0
        return ones(Int, n_ones)
    end

    core_format = config.core_format
    use_dense = core_format === :dense ||
                (core_format === :auto && nnz >= config.dense_min_density * nr * nc &&
                 8 * nr * nc <= config.dense_max_bytes)
    verbose && println("Core: nnz = $nnz (density $(round(100 * nnz / (nr * nc); digits = 2)) %), format ",
            use_dense ? "dense (ZZMatrix)" : "sparse (SMat)")

    t0 = time_ns()
    core = use_dense ? extract_core_dense(M) : extract_core_smat(M)
    divs_core = BigInt.(Hecke.elementary_divisors(core))
    if stats !== nothing
        stats.core_seconds = (time_ns() - t0) / 1e9
        stats.core_format = use_dense ? :dense : :sparse
    end

    divs_final = vcat(ones(Int, n_ones), divs_core)
    filter!(x -> x != 0, divs_final)
    return sort!(divs_final)
end

# The same from the rows, with the coefficients in value_type.
function _reduce_and_divisors(rows_data::RowsData{<:Integer}, m::Int, n::Int;
                              value_type::Type = Int32, kwargs...)
    M = DualSparseMatrix(rows_data, m, n; value_type)
    return _reduce_and_divisors!(M; kwargs...)
end

# A failed attempt of the widening, recorded in stats.
_record_failure!(::Nothing, T, e, seconds) = nothing
function _record_failure!(stats::ReductionStats, T, e::CoefficientOverflowError, seconds::Float64)
    push!(stats.failed_attempts, (value_type = T, unit_pivots = e.unit_pivots, seconds = seconds))
    return nothing
end

"""
    _reduce_and_divisors_widening(rows_data, m, n, value_types; config = ReductionConfig(),
                                  trace = nothing, stats = nothing, probe = nothing,
                                  verbose = false)

`_reduce_and_divisors` with each type of `value_types` in turn, from the input, until one
does not overflow (steps 2a and 2b of s3, v0.1.1). The trace and the counters of `stats`
are reset at the start of each attempt, so that they describe the attempt that succeeds
(the choice of the pivots does not depend on the type: the trajectory is that of a direct
computation with the type retained); `stats` also receives that type and the failed
attempts, the last one included when its overflow is rethrown. The retry happens outside
the `catch` block, so that nothing of a failed attempt stays reachable.
"""
function _reduce_and_divisors_widening(rows_data::RowsData{<:Integer}, m::Int, n::Int, value_types;
                                       config::ReductionConfig = ReductionConfig(),
                                       trace::Union{Nothing, ReductionTrace} = nothing,
                                       stats::Union{Nothing, ReductionStats} = nothing,
                                       probe = nothing, verbose::Bool = false)
    stats === nothing || _reset!(stats)
    k = length(value_types)
    for (i, T) ∈ enumerate(value_types)
        trace === nothing || _reset!(trace)
        stats === nothing || _reset_counters!(stats)
        t0 = time_ns()
        result = try
            _reduce_and_divisors(rows_data, m, n; value_type = T, config, trace, stats, probe, verbose)
        catch e
            e isa CoefficientOverflowError || rethrow()
            if i == k
                _record_failure!(stats, T, e, (time_ns() - t0) / 1e9)
                rethrow()
            end
            e
        end
        if result isa CoefficientOverflowError
            seconds = (time_ns() - t0) / 1e9
            _record_failure!(stats, T, result, seconds)
            verbose && println("value_type = $T: overflow after $(result.unit_pivots) unit pivots ",
                               "($(round(seconds; digits = 2)) s), retrying with $(value_types[i + 1])")
            continue
        end
        stats === nothing || (stats.value_type = T)
        verbose && println("value_type = $T retained")
        return result
    end
end

#############################
### Dense core extraction ###
#############################

# Motivation (ZT2w4m5 d₄): core 9351×25417, nnz = 130 M; extraction to a Hecke `SMat` took 422 s (as long as the whole reduction), i.e. ~3.2 µs/entry versus ~0.25 µs/entry on G₃₇ d₅ (nnz = 2.7 M). Probable cause (not proven): GC pressure (one `ZZRingElem` and one intermediate tuple allocated per entry, pointer vectors to traverse). The dense extraction below avoids these costs; it does not change the result of the reduction.

# A core coefficient for Flint: an Int when it fits a machine word (no allocation; for a
# ZZRingElem, Hecke.ZZ(::Int) avoids the BigInt that Hecke.ZZ(::Int32) goes through), a
# ZZRingElem otherwise. Int(v) on an Int128 beyond Int64 was the InexactError of v0.1.0.
@inline _to_zz(v::Union{Int8, Int16, Int32, Int64}) = Int(v)
@inline _to_zz(v::Int128) = typemin(Int) <= v <= typemax(Int) ? Int(v) : Hecke.ZZ(v)
@inline _to_zz(v::Integer) = Hecke.ZZ(v)

function _core_index_maps(M::DualSparseMatrix)
    active_rows = [r for r ∈ 1:M.num_rows if !isempty(M.row_cols[r])]
    active_cols = [c for c ∈ 1:M.num_cols if !isempty(M.col_rows[c])]
    col_inverse = zeros(Int, M.num_cols)
    for (new_idx, old_idx) ∈ enumerate(active_cols)
        col_inverse[old_idx] = new_idx
    end
    return active_rows, active_cols, col_inverse
end

"""
    extract_core_dense(M) -> ZZMatrix

The core (active rows × active columns, same numbering) as a DENSE Flint matrix: a single C allocation, no per-entry `ZZRingElem` object on the Julia GC side. Memory cost: 8 bytes × rows × columns (ZT2w4m5 d₄: 9351×25417 ≈ 1.9 GiB if the entries fit in a machine word).
"""
function extract_core_dense(M::DualSparseMatrix)
    active_rows, active_cols, col_inverse = _core_index_maps(M)
    nr, nc = length(active_rows), length(active_cols)
    A = Hecke.zero_matrix(Hecke.ZZ, nr, nc)
    (nr == 0 || nc == 0) && return A
    for (i, old_r) ∈ enumerate(active_rows)
        cols = M.row_cols[old_r]
        vals = M.row_vals[old_r]
        @inbounds for idx ∈ eachindex(cols)
            new_c = col_inverse[cols[idx]]
            new_c > 0 && (A[i, new_c] = _to_zz(vals[idx]))
        end
    end
    return A
end

