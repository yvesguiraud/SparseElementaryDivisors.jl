# Guide

```@meta
CurrentModule = SparseElementaryDivisors
```

This page explains how to call [`sparse_elementary_divisors`](@ref) and what its options do; the complete descriptions are in the [Reference](reference.md).

## Input

The matrix can be given as:
- a `SparseMatrixCSC` with integer entries (the fastest path), or any other `AbstractSparseMatrix` or `AbstractMatrix{<:Integer}`, which is converted first (dense matrices, adjoints and views work, but are not efficient);
- a Hecke matrix, `Hecke.SMat{ZZRingElem}` or `Hecke.ZZMatrix`;
- rows given one by one, `rows[r] = [(column, value), ...]` with the columns strictly increasing, together with the size of the matrix: `sparse_elementary_divisors(rows, m, n)` or `sparse_elementary_divisors(m, n, rows)` (see [`RowsData`](@ref SparseElementaryDivisors.RowsData));
- the name of a file in Matrix Market format (see [Matrix Market files](@ref matrix-market) below).

Entries of any integer type are accepted, including entries beyond a machine word (`BigInt`, `ZZRingElem`). Rows given one by one are checked before the computation: an `ArgumentError` names the first fault (column out of range, unsorted or repeated column, wrong number of rows). Unsorted triplets `(i, j, v)` can go through `sparse(I, J, V)` of SparseArrays, which sorts them and sums the duplicates.

## Result

```julia
using SparseElementaryDivisors, SparseArrays

A = sparse([2 0 0; 0 3 0; 0 0 0])
d = sparse_elementary_divisors(A)    # BigInt[1, 6]
```

The result is the vector of the non-zero elementary divisors, in increasing order, each one dividing the next: here `A` is equivalent to the diagonal matrix with entries 1, 6 and 0. Its length is the rank of the matrix; the remaining `min(m, n) - length(d)` elementary divisors are zero.

In a homology computation, with differentials ``d_n : C_{n+1} \to C_n``, the homology in degree ``n`` is ``\mathbb{Z}^{b} \oplus \bigoplus_i \mathbb{Z}/e_i\mathbb{Z}``, where ``b = \operatorname{rank} C_n - \operatorname{rank} d_{n-1} - \operatorname{rank} d_n`` and the ``e_i`` are the elementary divisors of ``d_n`` greater than 1.

## Threads

The pivot search and the application of the pivots run in parallel when Julia is started with several threads: `julia -t auto`, or `julia -t 16`. The keyword `nworkers` sets the number of parallel tasks, the calling thread included (by default `Threads.nthreads()`); it can be lowered on a shared machine, or for small matrices, where the parallel steps do not pay off. The pivots are chosen in the same order whatever the number of workers, so that the result does not depend on it.

## Integer types

During the elimination, the coefficients are machine integers with checked arithmetic. With the default `value_type = :auto`, the types of [`DEFAULT_VALUE_TYPES`](@ref SparseElementaryDivisors.DEFAULT_VALUE_TYPES), `(Int32, Int64, Int128, BigInt)`, are tried in turn, starting with the first one that holds all the entries; after an overflow, the computation restarts from the input with the next type. Each call is independent: a failed attempt costs its own time, and the next matrix starts again from the first type that holds its entries. Wider types cost more time and memory, `BigInt` much more, hence the order.

```julia
st = ReductionStats()
sparse_elementary_divisors(A; stats = st)
st.value_type          # the type retained
st.failed_attempts     # the attempts that overflowed: (type, unit pivots, seconds)

sparse_elementary_divisors(A; value_type = Int64)             # Int64 only
sparse_elementary_divisors(A; value_type = (Int64, Int128))   # Int64, then Int128
```

A fixed type or a sequence restricts the choice; a [`CoefficientOverflowError`](@ref) is raised when the last type does not suffice. Signed integers of fixed width beyond `Int128`, such as `Int256` of [BitIntegers.jl](https://github.com/rfourquet/BitIntegers.jl), are accepted as well.

## Settings

The settings of the reduction are gathered in a [`ReductionConfig`](@ref), an immutable value that can be kept and reused, for instance for all the differentials of a complex. Each of its fields can also be given directly as a keyword of `sparse_elementary_divisors`, and then takes precedence over `config`; an unknown keyword raises an `ArgumentError`.

```julia
cfg = ReductionConfig(nworkers = 4, core_format = :sparse)
sparse_elementary_divisors(A; config = cfg)
sparse_elementary_divisors(A; config = cfg, nworkers = 1)     # nworkers = 1 wins
```

No setting changes the elementary divisors. The parameters of the pivot search change the order of the pivots, hence the time, the fill-in and the size of the coefficients; the other settings only change how the work is done. The default values come from measurements on matrices like those of the [Benchmarks](benchmarks.md). Besides `nworkers`, the setting most often useful is `core_format`, the format of the core passed to Hecke: `:sparse`, `:dense`, or `:auto`, which chooses from the density and the size of the core.

## Statistics

A [`ReductionStats`](@ref) passed with the keyword `stats` is filled in place: the integer type retained and the attempts that overflowed, the numbers of unit pivots and of pivot searches, the times of the phases (structured elimination, pivot search, pivot application, core), and the size and format of the core. With `verbose = true`, the attempts that overflowed, the type retained and the size and format of the core are also printed.

```julia
st = ReductionStats()
sparse_elementary_divisors(A; stats = st, verbose = true)
st.unit_pivots, st.core_size, st.core_format
```

## [Matrix Market files](@id matrix-market)

`sparse_elementary_divisors(path)` reads the file `path` with [`read_matrix_market`](@ref), which accepts the [Matrix Market](https://math.nist.gov/MatrixMarket/formats.html) formats `coordinate` and `array`, of type `integer` or `pattern` (and `real` when all the values are integers), `general`, `symmetric` or `skew-symmetric`. [`write_matrix_market`](@ref) writes a matrix in the format `coordinate integer general`, entries sorted by row, then by column.

```julia
write_matrix_market("A.mtx", A; comment = "a small example")
B = read_matrix_market("A.mtx")
```

The values are read as `Int`: an entry beyond a machine word can be written, but not read back yet.

## Comparison with Magma

When Magma is available, locally or on a machine reachable by `ssh`, [`magma_elementary_divisors`](@ref) computes the elementary divisors of the same matrix with Magma's `ElementaryDivisors`, and [`timed_magma_elementary_divisors`](@ref) also returns the time that Magma measures for that call alone. These functions serve to check results and compare times; Magma is not needed otherwise.

```julia
divs = magma_elementary_divisors(A; ssh_host = "myserver")      # or locally: magma_elementary_divisors(A)
divs, t = timed_magma_elementary_divisors(A; ssh_host = "myserver")
```

## Method

The matrix is stored both by rows and by columns, and reduced in three steps:
1. a structured elimination (LaMacchia and Odlyzko) pivots on the rows and columns that have a single entry, when this entry is ±1;
2. each following step chooses an entry ±1 of small Markowitz score ``(r - 1)(c - 1)``, where ``r`` and ``c`` are the numbers of entries of its row and of its column, eliminates its row and its column, and runs the structured elimination again;
3. when no entry ±1 is left, the remaining matrix, the core, goes to `Hecke.elementary_divisors` (Hecke and Nemo, on FLINT), as a sparse or a dense matrix.

Each pivot ±1 gives an elementary divisor 1. The method suits matrices with many entries ±1 and few entries per row and column, such as the differentials of resolutions, where the core is small. When few pivots ±1 are available, most of the work is done by Hecke on a large core. The articles behind these steps are listed in [References](references.md).
