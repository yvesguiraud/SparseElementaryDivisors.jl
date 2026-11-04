# References

## Methods

- H. M. Markowitz, *The Elimination Form of the Inverse and its Application to Linear Programming*, Management Science 3(3) (1957), 255–269, [doi:10.1287/mnsc.3.3.255](https://doi.org/10.1287/mnsc.3.3.255). The Markowitz score (r − 1)(c − 1) that orders the search of the unit pivots.
- B. A. LaMacchia, A. M. Odlyzko, *Solving Large Sparse Linear Systems Over Finite Fields*, in A. J. Menezes, S. A. Vanstone (eds.), Advances in Cryptology – CRYPTO '90, Lecture Notes in Computer Science 537, Springer (1991), 109–133, [doi:10.1007/3-540-38424-3_8](https://doi.org/10.1007/3-540-38424-3_8). Structured elimination, here of the rows and columns of degree one with a unit entry, before and between the Markowitz pivots.
- J.-G. Dumas, B. D. Saunders, G. Villard, *On Efficient Sparse Integer Matrix Smith Normal Form Computations*, Journal of Symbolic Computation 32(1–2) (2001), 71–99, [doi:10.1006/jsco.2001.0451](https://doi.org/10.1006/jsco.2001.0451). Smith normal form of large sparse integer matrices, in particular from homology.
- J.-G. Dumas, F. Heckenbach, B. D. Saunders, V. Welker, *Computing Simplicial Homology Based on Efficient Smith Normal Form Algorithms*, in M. Joswig, N. Takayama (eds.), Algebra, Geometry and Software Systems, Springer (2003), 177–206, [doi:10.1007/978-3-662-05148-1_10](https://doi.org/10.1007/978-3-662-05148-1_10). Review and comparison of Smith normal form algorithms on sparse boundary matrices.

## Software used

- C. Fieker, W. Hart, T. Hofmann, F. Johansson, *Nemo/Hecke: Computer Algebra and Number Theory Packages for the Julia Programming Language*, in Proceedings of ISSAC '17, ACM (2017), 157–164, [doi:10.1145/3087604.3087611](https://doi.org/10.1145/3087604.3087611). Hecke computes the elementary divisors of the core.
- The FLINT team, *FLINT: Fast Library for Number Theory* (2025), [flintlib.org](https://flintlib.org). The integer arithmetic of the core, through Nemo.

## Comparison with Magma

- W. Bosma, J. Cannon, C. Playoust, *The Magma Algebra System I: The User Language*, Journal of Symbolic Computation 24(3–4) (1997), 235–265, [doi:10.1006/jsco.1996.0125](https://doi.org/10.1006/jsco.1996.0125).

## File format

- R. F. Boisvert, R. Pozo, K. A. Remington, *The Matrix Market Exchange Formats: Initial Design*, NIST IR 5935, National Institute of Standards and Technology (1996), [doi:10.6028/NIST.IR.5935](https://doi.org/10.6028/NIST.IR.5935). The format of [`read_matrix_market`](@ref SparseElementaryDivisors.read_matrix_market) and [`write_matrix_market`](@ref SparseElementaryDivisors.write_matrix_market).

## Matrices of the benchmarks

- P. Dehornoy, Y. Lafont, *Homology of Gaussian Groups*, Annales de l'Institut Fourier 53(2) (2003), 489–540, [doi:10.5802/aif.1951](https://doi.org/10.5802/aif.1951). The resolution whose reduced version gives the G₃₇ matrices.
- N. Idrissi, V. Roca i Lucio, *Homology of configuration spaces in positive characteristic via point-set constructions*, [arXiv:2606.26802](https://arxiv.org/abs/2606.26802). The ZT2w4m5 matrices.
