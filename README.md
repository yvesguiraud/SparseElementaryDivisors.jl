# SparseElementaryDivisors.jl

[![CI](https://github.com/yvesguiraud/SparseElementaryDivisors.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yvesguiraud/SparseElementaryDivisors.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/yvesguiraud/SparseElementaryDivisors.jl/graph/badge.svg)](https://codecov.io/gh/yvesguiraud/SparseElementaryDivisors.jl)

Computation of the elementary divisors of large sparse integer matrices, in Julia.

The matrix is first reduced by exact sparse elimination with unit pivots (Markowitz-style pivot choice, parallel pivot search and parallel fill-in application, checked integer arithmetic); then the small remaining matrix, the *core*, is handed to [Hecke.jl](https://github.com/thofma/Hecke.jl). The computation is multithreaded, but its result does not depend on the number of threads used.

The package was initially developed for the differentials of resolutions of Garside monoids (mostly in the context of the homology of real and complex braid groups), whose matrices have hundreds of thousands of rows and columns. One aim is to bring to Julia performances similar to those of [Magma](https://magma-maths.org) on such matrices.

> **Status:** early release (0.1); the interface may still change.

The [documentation](https://yvesguiraud.github.io/SparseElementaryDivisors.jl/stable/) contains a user guide, the benchmarks, the reference manual and a bibliography, as well as the development history, the acknowledgements and how to cite the package.

## Installation

```julia
using Pkg
Pkg.add("SparseElementaryDivisors")
```

or press `]` to enter the package mode of the REPL, then `add SparseElementaryDivisors`.

## First examples

```julia
using SparseElementaryDivisors, SparseArrays

A = sparse([2 0 0; 0 3 0; 0 0 0])
sparse_elementary_divisors(A)              # returns BigInt[1, 6]

sparse_elementary_divisors("matrix.mtx")   # a matrix in Matrix Market format
```

Start Julia with several threads to benefit from the parallel steps: `julia -t auto`. The [user guide](https://yvesguiraud.github.io/SparseElementaryDivisors.jl/stable/user_guide/) describes the accepted inputs and the options.

## Benchmarks

The [benchmarks](https://yvesguiraud.github.io/SparseElementaryDivisors.jl/stable/benchmarks/) compare the times of the package with those of Hecke.jl and Magma on large matrices from homology computations.
