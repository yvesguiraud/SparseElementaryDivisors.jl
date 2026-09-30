# SparseElementaryDivisors.jl

Elementary divisors (Smith normal form) of **large sparse integer matrices**, in Julia.

The matrix is first reduced by exact sparse elimination with unit pivots (Markowitz-style pivot choice, parallel pivot search and parallel fill-in application, checked integer arithmetic); the small remaining "core" is handed to [Hecke.jl](https://github.com/thofma/Hecke.jl). The result does not depend on the number of threads. The package was developed for the differentials of resolutions of Garside monoids (mostly in the context of the homology of real and complex braid groups), whose matrices have hundreds of thousands of rows and columns, and aims at the performance of Magma on such inputs.

Any positive or negative experience on other matrices is most welcome!

> **Status: early release (0.1).** The API may still change.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/yvesguiraud/SparseElementaryDivisors.jl")
```

## Usage

```julia
using SparseElementaryDivisors, SparseArrays

A = sparse([2 0 0; 0 3 0; 0 0 0])
sparse_elementary_divisors(A)                 # BigInt[1, 6]  (the non-zero elementary divisors)

sparse_elementary_divisors("matrix.mtx")      # a Matrix Market file
```

Start Julia with several threads to benefit from the parallel steps: `julia -t auto`.

Accepted inputs:
- `SparseMatrixCSC` (fast path), any other `AbstractSparseMatrix` or `AbstractMatrix{<:Integer}` (converted; dense matrices, adjoints and views work but are not efficient);
- Hecke matrices: `Hecke.SMat{ZZRingElem}` and `Hecke.ZZMatrix`;
- the matrix given row by row, `rows[r] = [(column, value), ...]` (sorted by column), together with its size: `sparse_elementary_divisors(rows, m, n)` or `sparse_elementary_divisors(m, n, rows)`;
- a file name: the matrix is read in the [Matrix Market](https://math.nist.gov/MatrixMarket/formats.html) format (`coordinate` or `array`, `integer` or `pattern`, `general`, `symmetric` or `skew-symmetric`) by `read_matrix_market`, which can also be called directly.

Main options:

| keyword | meaning |
|---|---|
| `value_type = Int32` | integer type used during the elimination; an arithmetic overflow raises an explicit error, in which case retry with `Int64`, `Int128` or `BigInt` |
| `nworkers` | number of parallel tasks (default: `Threads.nthreads()`; lower it if the machine is shared, or for tiny matrices) |
| `core_format = :auto` | format of the core passed to Hecke: `:sparse`, `:dense` or `:auto` |
| `verbose = false` | print the size and format of the core |

## Comparison with Magma

If the user has access to Magma (locally, or through `ssh`), `magma_elementary_divisors` runs Magma's `ElementaryDivisors` on the same matrix and returns the divisors. The variant `timed_magma_elementary_divisors` also returns Magma's own computation time:

```julia
divs = magma_elementary_divisors(A; ssh_host = "myserver")          # or locally: magma_elementary_divisors(A)
divs, t = timed_magma_elementary_divisors(A; ssh_host = "myserver")
```

## Benchmarks (informational)

Times for the complete computation of the elementary divisors (for SparseElementaryDivisors: reduction, extraction of the core, and Hecke on the core), on **three different machines**:
- plain Hecke: 13th Gen Intel Core i9-13900, 128 GB;
- SparseElementaryDivisors: Intel Core Ultra 9 285 vPro, 128 GB;
- Magma (with Magma's own timing of `ElementaryDivisors`): Intel Xeon Gold 5222, 64 GB.

The times are only indicative: the test matrices are not distributed with the package (this may change), and the machines are not exactly comparable. The matrices are differentials of a reduced version of Dehornoy–Lafont's resolution for the dual braid monoid of the complex reflection group G₃₇, produced by [Gauss.jl](https://plmlab.math.cnrs.fr/guiraud/gauss.jl), and a family ZT2w4m5 of examples provided by Najib Idrissi (see Acknowledgements).

| Matrix | Size | Hecke | SparseElementaryDivisors | Magma |
|---|---|---|---|---|
| G₃₇ d₃ | 15120 × 54327 | 4 min | 2.2 s | 5.5 s |
| G₃₇ d₄ | 54327 × 108360 | 10 h 49 min | 11.3 s | ≈ 1 min |
| G₃₇ d₅ | 108360 × 121555 | 13 days 21 h | 23.1 s | ≈ 45 s |
| G₃₇ d₆ | 121555 × 71760 | — | 12.2 s | ≈ 20 s |
| ZT2w4m5 d₂ | 43384 × 166406 | — | 2.7 s | 6.3 s |
| ZT2w4m5 d₃ | 166406 × 562046 | — | 24.0 s | 53 s |
| ZT2w4m5 d₄ | 562046 × 1776117 | — | 6 min 34 s | 24 min 30 s |

Notes: the SparseElementaryDivisors times are medians of 3 to 5 runs (2 runs for ZT2w4m5 d₄) of the whole computation, using all 24 threads (`nworkers = Threads.nthreads()`). The final step on the core by Hecke takes less than 0.4 s in all cases except ZT2w4m5 d₄, where it takes about 19 s (the core has size 9351 × 25417 and is dense).

## Development

The package was developed with the help of an AI assistant.

The starting point was a standard request to Google in August 2026, to see which tools existed to compute the elementary divisors of huge sparse matrices. Google's AI mode had just appeared and, starting from its answer, the conversation progressively led to a plan to improve the capabilities of Hecke.jl itself. It produced a first version that was already an improvement: for example, it reduced the computation for d₅ in type G₃₇ from almost two weeks to 45 minutes. At that point, Google's AI mode indicated that it was using the Claude Sonnet 3.5 model.

After that first prototype, the author used Claude (with the Sonnet 5, then Sonnet 5.5 models) to help with the subsequent development of the code, with the design of the documentation, and with the finalisation of the package.

## Acknowledgements

The author thanks Najib Idrissi for providing the family of examples called ZT2w4m5 here, which are resolutions linked to the homology of configuration spaces on the torus. These examples greatly helped the development of this package. They come from the following article: Najib Idrissi and Victor Roca i Lucio, *Homology of configuration spaces in positive characteristic via point-set constructions*, [arXiv:2606.26802](https://arxiv.org/abs/2606.26802).

## License

MIT.
