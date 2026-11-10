# Benchmarks

Times for the complete computation of the elementary divisors (for SparseElementaryDivisors: reduction, extraction of the core, and Hecke on the core), on **three different machines**:
- plain Hecke: 13th Gen Intel Core i9-13900, 128 GB;
- SparseElementaryDivisors: Intel Core Ultra 9 285 vPro, 128 GB;
- Magma (with Magma's own timing of `ElementaryDivisors`): Intel Xeon Gold 5222, 64 GB.

The times are only indicative: the test matrices are not distributed with the package (this may change), and the machines are not exactly comparable. The matrices are differentials of a reduced version of Dehornoy–Lafont's resolution (see the [Bibliography](bibliography.md)) for the dual braid monoid of the complex reflection group G₃₇, produced by the preliminary implementation of [GarsideHomology.jl](https://github.com/yvesguiraud/GarsideHomology.jl), and a family ZT2w4m5 of examples provided by Najib Idrissi (see the [acknowledgements](@ref acknowledgements)).

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
