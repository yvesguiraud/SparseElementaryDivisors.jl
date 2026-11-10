# SparseElementaryDivisors.jl

Computation of the elementary divisors of large sparse integer matrices, in Julia.

!!! note "Status"
    Early release (0.1): the interface may still change.

## Introduction

### Elementary divisors and homology

For every matrix ``A`` over a principal ideal domain ``R``, there exist invertible matrices ``S`` and ``T`` over ``R`` such that

```math
SAT = \begin{pmatrix} D & 0 \\ 0 & 0 \end{pmatrix},
\qquad D = \begin{pmatrix} d_1 & & 0 \\ & \ddots & \\ 0 & & d_r \end{pmatrix},
```

where ``d_1, \dots, d_r`` are non-zero elements of ``R``, and each ``d_i`` divides ``d_{i+1}``; moreover, these diagonal elements are unique up to multiplication by units of ``R``. The diagonal matrix ``SAT`` is called the *Smith normal form* of ``A``, and ``d_1, \dots, d_r`` are called the *elementary divisors* (or *invariant factors*) of ``A``.

A common use of elementary divisors is the computation of homology groups. For example, consider a complex

```math
\mathbb{Z}^m \xrightarrow{f} \mathbb{Z}^n \xrightarrow{g} \mathbb{Z}^p,
```

that is, with ``gf = 0``. Let ``r`` and ``s`` be the respective ranks of ``f`` and ``g``, and ``d_1, \dots, d_r`` be the elementary divisors of ``f``. Then the homology of the complex at the middle term is

```math
\ker g / \operatorname{im} f \;\cong\; \mathbb{Z}^{\,n - r - s} \oplus \bigoplus_{i=1}^{r} \mathbb{Z}/d_i\mathbb{Z}.
```

### Computation method

The package computes *positive* elementary divisors of integer matrices, which makes the result unique. The method is specially designed for *sparse* matrices.

The matrix is first reduced by exact sparse elimination with unit pivots (Markowitz-style pivot choice, parallel pivot search and parallel fill-in application, checked integer arithmetic). Then the small remaining matrix, the *core*, is handed to [Hecke.jl](https://github.com/thofma/Hecke.jl). The computation is multithreaded, but its result does not depend on the number of threads used. The [User guide](user_guide.md) describes the method in more detail.

### Motivation

The package was initially developed for the differentials of resolutions of Garside monoids (mostly in the context of the homology of real and complex braid groups), whose matrices have hundreds of thousands of rows and columns. One aim is to bring to Julia performances similar to those of [Magma](https://magma-maths.org) on such matrices; see [Benchmarks](benchmarks.md) for current performance.

## Quick start

### Installation

```julia
using Pkg
Pkg.add("SparseElementaryDivisors")
```

or press `]` to enter the package mode of the REPL, then `add SparseElementaryDivisors`.

### First examples

The main function of the package is [`sparse_elementary_divisors`](@ref), which returns the elementary divisors of an integer matrix.

The argument can be given in various formats, such as a sparse integer matrix of type `SparseMatrixCSC`:

```julia
using SparseElementaryDivisors, SparseArrays

A = sparse([2 0 0; 0 3 0; 0 0 0])
sparse_elementary_divisors(A)            # returns BigInt[1, 6]
```

The argument can also be the name of a file that contains an integer matrix in Matrix Market format:

```julia
using SparseElementaryDivisors

sparse_elementary_divisors("matrix.mtx")
```

Start Julia with several threads to benefit from the parallel steps: `julia -t auto`. The [User guide](user_guide.md) describes the accepted inputs and the options; the [Reference manual](reference_manual.md) gathers the docstrings.

## About

### Development

This package is developed with the help of an AI assistant. In every version, the part of the code that is produced by the AI assistant is reviewed by the author and validated by tests against reference computations.

The starting point was a standard request to Google in August 2026, about the existing Julia tools to compute the elementary divisors of huge sparse matrices. Google's AI mode had just appeared and, starting from its answer, the conversation progressively led to a first prototype. In that first version, for example, the computation of the elementary divisors of the differential d₅ in type G₃₇ was reduced from almost two weeks with Hecke.jl to 45 minutes. At that point, Google's AI mode indicated that it was using the Claude Sonnet 3.5 model.

After the first prototype, the author used Claude (with the Sonnet 5, then Sonnet 5.5 models) to help with the subsequent development of the code, with the design of the documentation, and with the finalisation of version 0.1.0 of the package. The development continues, now using the Claude Opus 5.5 model.

### Citing

To cite SparseElementaryDivisors.jl, see the file [CITATION.cff](https://github.com/yvesguiraud/SparseElementaryDivisors.jl/blob/main/CITATION.cff) of the repository, or "Cite this repository" on its GitHub page. The methods and the software it relies on are listed in the [Bibliography](bibliography.md).

### [Acknowledgements](@id acknowledgements)

The author thanks Jean Michel for having introduced him to Julia a few years ago, and for his regular and valuable advice on the language. His package [Chevie.jl](https://github.com/jmichel7/Chevie.jl) is particularly useful for generating concrete test examples.

The author also thanks Najib Idrissi for providing the family of examples called ZT2w4m5 here, which are resolutions linked to the homology of configuration spaces on the torus (from the article: Najib Idrissi and Victor Roca i Lucio, *Homology of configuration spaces in positive characteristic via point-set constructions*, [arXiv:2606.26802](https://arxiv.org/abs/2606.26802)). These examples helped the development of this package.

### Licence

MIT, see the file [LICENSE](https://github.com/yvesguiraud/SparseElementaryDivisors.jl/blob/main/LICENSE) of the repository.

### Contact

Feedback on other matrices is most welcome, by [email](mailto:yves.guiraud@imj-prg.fr) to the author, or as an [issue](https://github.com/yvesguiraud/SparseElementaryDivisors.jl/issues) on the repository.
