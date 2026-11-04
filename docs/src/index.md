# SparseElementaryDivisors.jl

Elementary divisors (Smith normal form) of **large sparse integer matrices**, in Julia.

The matrix is first reduced by exact sparse elimination with unit pivots (Markowitz-style pivot choice, parallel pivot search and parallel fill-in application, checked integer arithmetic); the small remaining "core" is handed to [Hecke.jl](https://github.com/thofma/Hecke.jl). The result does not depend on the number of threads. The package was developed for the differentials of resolutions of Garside monoids (mostly in the context of the homology of real and complex braid groups), whose matrices have hundreds of thousands of rows and columns, and aims at the performance of Magma on such inputs.

Any positive or negative experience on other matrices is most welcome, for instance as an [issue](https://github.com/yvesguiraud/SparseElementaryDivisors.jl/issues) of the repository.

!!! note "Status"
    Early release (0.1): the interface may still change.

## Installation

```julia
using Pkg
Pkg.add("SparseElementaryDivisors")
```

or press `]` to enter the package mode of the REPL, then `add SparseElementaryDivisors`.

## A first example

```julia
using SparseElementaryDivisors, SparseArrays

A = sparse([2 0 0; 0 3 0; 0 0 0])
sparse_elementary_divisors(A)            # BigInt[1, 6], the non-zero elementary divisors

sparse_elementary_divisors("matrix.mtx") # a file in Matrix Market format
```

Start Julia with several threads to benefit from the parallel steps: `julia -t auto`. The [Guide](guide.md) describes the accepted inputs and the options; the [Reference](reference.md) gathers the docstrings.

## Citing

To cite SparseElementaryDivisors.jl, see the file [CITATION.cff](https://github.com/yvesguiraud/SparseElementaryDivisors.jl/blob/main/CITATION.cff) of the repository, or "Cite this repository" on its GitHub page. The methods and the software it relies on are listed in [References](references.md).

## Development

The package was developed with the help of an AI assistant.

The starting point was a standard request to Google in August 2026, to see which tools existed to compute the elementary divisors of huge sparse matrices. Google's AI mode had just appeared and, starting from its answer, the conversation progressively led to a plan to improve the tools to compute elementary divisors in Julia: starting from the functions of Hecke.jl, the goal is to try to reach performances similar to those of Magma, but in an open world. A first version was produced, and it was already an improvement: for example, it reduced the computation for d₅ in type G₃₇ from almost two weeks (with Hecke.jl) to 45 minutes. At that point, Google's AI mode indicated that it was using the Claude Sonnet 3.5 model.

After that first prototype, the author used Claude (with the Sonnet 5, then Sonnet 5.5 models) to help with the subsequent development of the code, with the design of the documentation, and with the finalisation of version 0.1.0 of the package. The development continues, now using the Claude Opus 5.5 model.

## [Acknowledgements](@id acknowledgements)

The author thanks Najib Idrissi for providing the family of examples called ZT2w4m5 here, which are resolutions linked to the homology of configuration spaces on the torus. These examples greatly helped the development of this package. They come from the following article: Najib Idrissi and Victor Roca i Lucio, *Homology of configuration spaces in positive characteristic via point-set constructions*, [arXiv:2606.26802](https://arxiv.org/abs/2606.26802).

## License

MIT, see the file [LICENSE](https://github.com/yvesguiraud/SparseElementaryDivisors.jl/blob/main/LICENSE) of the repository.
