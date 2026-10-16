# Test data

Matrices in Matrix Market format (`coordinate integer general`) and reference results
used by `test/test_fingerprints.jl`.

## Matrices

- `G34_d<i>.mtx`, `G35_d<i>.mtx`, `G36_d<i>.mtx`: the differentials dᵢ : Cᵢ₊₁ → Cᵢ
  (rows indexed by Cᵢ) of the integral chain complexes computed by Gauss.jl for the
  dual braid monoids of the complex reflection groups G₃₄, G₃₅ and G₃₆.
- `synth_*.mtx`: synthetic matrices generated once with a fixed seed (the comment line
  of each file describes it). `synth_b_150x200` overflows `Int32` during elimination;
  `synth_snf_100x120` has known elementary divisors; `synth_dense_60x60` has no
  acceptable unit pivot and goes to the dense core almost at once.

The differentials of G₃₇ (about 5.8 million nonzero entries) are not stored here. Their
references are in `references.toml` with `heavy = true`; the corresponding tests run only
when the environment variable `SED_HEAVY_DIR` gives the folder of their `.mtx` files.

## References

`references.toml` was computed with SparseElementaryDivisors v0.1.0 (Julia 1.13.1; the
same values were obtained with Julia 1.13.2). For each matrix: rank, number of unit
divisors, torsion (as `divisor:multiplicity`), value type used (`Int32`, else `Int64`,
else `BigInt`, with the list of types that overflowed), trajectory fingerprint of the
reduction (`trace_hash`, `trace_pivots`, `trace_sge` for the number of structured elimination
passes, `trace_sge_units` for the units they found) and size of the core. The fingerprints were checked to be
identical with the default configuration and with the parallel code paths forced, for
1, 2 and 4 workers.
