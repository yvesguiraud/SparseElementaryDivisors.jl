# Benchmarks

Times of the complete computation of the elementary divisors of the matrices below, on three machines. For the package, the time is that of the whole call to [`sparse_elementary_divisors`](@ref) on the matrix already read from its file: reduction, extraction of the core, and Hecke on the core. For Magma, it is the time that Magma measures for `ElementaryDivisors` alone, without the transfer of the matrix. For Hecke alone, it is the time of `Hecke.elementary_divisors` on the matrix already built as a `SMat`. Each time is the median of three runs, unless stated otherwise. The measurements were made in November 2026 with version 0.1.1 of the package, Julia 1.13.3, Hecke 0.41.1 and Magma V2.29-4, by the script `benchmark/benchmarks.jl` of the repository, where the raw results are kept (`benchmark/results`).

## Matrices

The matrices G₃₇ d₃ to d₆ are the differentials of the reduced lcm resolution of the dual braid monoid of the complex reflection group G₃₇ (a variant of Dehornoy–Lafont's resolution with fewer cells, see the [Bibliography](bibliography.md)), produced by the preliminary implementation of [GarsideHomology.jl](https://github.com/yvesguiraud/GarsideHomology.jl). The matrices ZT2w4m5 d₂ to d₄ were provided by Najib Idrissi (see the [acknowledgements](@ref acknowledgements)). They are not distributed with the package for now. The core is the matrix left after the elimination of the unit pivots, which Hecke then reduces.

| Matrix | Size | Entries | Core | Elementary divisors |
|---|---|---|---|---|
| G₃₇ d₃ | 15120 × 54327 | 290302 | 1 × 30790 | 1 (13103 times), 2 |
| G₃₇ d₄ | 54327 × 108360 | 967056 | 15 × 57749 | 1 (41221 times), 2, 6 |
| G₃₇ d₅ | 108360 × 121555 | 1778659 | 253 × 40710 | 1 (67135 times), 3, 6 |
| G₃₇ d₆ | 121555 × 71760 | 1830531 | 529 × 6018 | 1 (54416 times), 2, 6 |
| ZT2w4m5 d₂ | 43384 × 166406 | 1235520 | 12 × 3713 | 1 (35526 times), 2 |
| ZT2w4m5 d₃ | 166406 × 562046 | 4883524 | 2494 × 14180 | 1 (130873 times), 2 (3 times) |
| ZT2w4m5 d₄ | 562046 × 1776117 | 17367804 | 9351 × 25417 | 1 (431162 times), 2 (3 times) |

## Machines

| | Processor | Threads | Memory | System |
|---|---|---|---|---|
| Server 1 | Intel Core i9-13900 | 32 | 128 GB | Linux |
| Server 2 | 2 × Intel Xeon Gold 5222 | 16 | 64 GB | Linux |
| Workstation | Intel Core Ultra 9 285H | 16 | 64 GB | Windows |

Server 1 is shared: jobs of other users were running during the measurements; the runs of a same series differed by at most 11%. On server 2, a computation on one thread takes about 25% longer on one of its two processors than on the other, and the system runs it on either; its times on one thread vary accordingly. Magma is available to us on server 2 only, so the comparison with Magma is made there, on the same matrices and the same machine. The package runs wherever Julia runs: the other tables give its times on server 1, with Hecke alone for comparison, and on a workstation.

## Comparison with Hecke

On server 1: Hecke alone, run on the first three matrices only, and the package with one thread and with all 32 threads.

| Matrix | Hecke alone | Package, 1 thread | Package, 32 threads |
|---|---|---|---|
| G₃₇ d₃ | 2 min 56 s | 3.2 s | 2.2 s |
| G₃₇ d₄ | 9 h 3 min | 20.5 s | 12.0 s |
| G₃₇ d₅ | 13 days 21 h [^1] | 1 min 23 s | 36.5 s |
| G₃₇ d₆ | — | 39.1 s | 16.4 s |
| ZT2w4m5 d₂ | — | 9.3 s | 3.8 s |
| ZT2w4m5 d₃ | — | 1 min 50 s | 24.2 s |
| ZT2w4m5 d₄ | — | 28 min 57 s [^2] | 4 min 54 s [^2] |

[^1]: Measured in August 2026 with Hecke 0.39.22, on the same machine.

[^2]: Median of two runs only.

## Comparison with Magma

On server 2: Magma, which uses one thread by default, and the package with one thread and with all 16 threads.

| Matrix | Magma | Package, 1 thread | Package, 16 threads |
|---|---|---|---|
| G₃₇ d₃ | 5.2 s | 6.4 s | 4.5 s |
| G₃₇ d₄ | 42.2 s | 46.9 s | 27.5 s |
| G₃₇ d₅ | 35.5 s | 2 min 58 s | 1 min 29 s |
| G₃₇ d₆ | 17.7 s | 1 min 27 s | 41.0 s |
| ZT2w4m5 d₂ | 5.5 s [^3] | 20.0 s | 7.2 s |
| ZT2w4m5 d₃ | 42.6 s [^3] | 3 min 46 s | 53.2 s |
| ZT2w4m5 d₄ | 28 min 45 s [^4] | 43 min 27 s [^2] | 11 min 31 s [^2] |

[^3]: Median of six runs.

[^4]: Median of six runs; see below.

Both computations have two phases. For a sparse matrix, Magma's `ElementaryDivisors` first performs a sparse elimination with Markowitz pivoting, then passes the remaining matrix to its algorithm for dense matrices, usually a modular one (see the chapter [Sparse Matrices](https://magma.maths.usyd.edu.au/magma/handbook/sparse_matrices) and the section [Canonical Forms](https://magma.maths.usyd.edu.au/magma/handbook/text/280) of the Magma handbook). The package eliminates the unit pivots, then passes the remaining matrix to Hecke. In both cases, the matrix passed to the second phase is called the core below.

Magma's `ElementaryDivisors` draws random numbers (the step of the random number generator, given by `GetSeed`, advances by millions during the call), from a seed that Magma chooses at random when it starts. Its core, and with it its time and its memory, change from one run to the next, as its verbose output shows (`SetVerbose("SparseMatrix", 1)`): on ZT2w4m5 d₄, the time of Magma ranges from 21 min 51 s to 56 min 17 s, and its peak memory, in the runs where it was measured, from 40 GiB to more than 59 GiB; two runs out of eight were stopped by the system for lack of memory. With a fixed seed (`magma -S 1`, see [Random Object Generation](https://magma.maths.usyd.edu.au/magma/handbook/text/18) in the Magma handbook), the computation is reproducible: on ZT2w4m5 d₄, 24 to 25 min and 26 GiB with this seed. The package draws no random numbers: its core is the same on every run, on every machine and for any number of threads. Its peak memory on ZT2w4m5 d₄ is 30 to 32 GiB on server 2.

### Size of the core

| Matrix | Core of Magma | Entries | Core of the package | Entries |
|---|---|---|---|---|
| G₃₇ d₃ | 168 × 1499 | 132332 | 1 × 30790 | 30790 |
| G₃₇ d₄ | 599 × 44702 | 1764035 | 15 × 57749 | 532291 |
| G₃₇ d₅ | 229 × 41866 | 1759603 | 253 × 40710 | 2709762 |
| G₃₇ d₆ | 417 × 2082 | 709618 | 529 × 6018 | 2564729 |
| ZT2w4m5 d₂ | 1267 × 3510 | 2238970 | 12 × 3713 | 44556 |
| ZT2w4m5 d₃ | 1856 × 21867 | 331682 | 2494 × 14180 | 240620 |
| ZT2w4m5 d₄ | 22342 × 28744 | 2765178 | 9351 × 25417 | 130357078 |

The core of Magma depends on the seed: the table gives it for `magma -S 1`, from the verbose output, in the orientation of the original matrix. The two eliminations choose different pivots and stop at different points. On ZT2w4m5 d₂, Magma stops when the density of the remaining matrix exceeds 50%, while the package continues. On ZT2w4m5 d₄, after almost the same number of pivots (431161 for Magma, 431162 for the package), the core of Magma has 2.8 million entries and that of the package 130 million.

## Times on a workstation

On the workstation: the package with one thread and with all 16 threads, both series running in the same process, and the peak memory of that process.

| Matrix | 1 thread | 16 threads | Peak memory |
|---|---|---|---|
| G₃₇ d₃ | 3.4 s | 2.4 s | 0.9 GiB |
| G₃₇ d₄ | 25.8 s | 14.0 s | 1.6 GiB |
| G₃₇ d₅ | 1 min 39 s | 26.6 s | 3.8 GiB |
| G₃₇ d₆ | 47.9 s | 13.9 s | 2.7 GiB |
| ZT2w4m5 d₂ | 10.4 s | 3.4 s | 1.6 GiB |
| ZT2w4m5 d₃ | 2 min 33 s | 33.0 s | 5.1 GiB |
| ZT2w4m5 d₄ | 50 min 21 s [^2] | 8 min 57 s [^2] | 33.3 GiB |

The final step, Hecke on the core, takes less than a second in all cases except ZT2w4m5 d₄, whose core is dense: 13 to 33 s depending on the machine.
