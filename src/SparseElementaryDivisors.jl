"""
    SparseElementaryDivisors

Elementary divisors (Smith normal form invariants) of large **sparse integer matrices**, by exact sparse elimination (unit pivots with Markowitz-style ordering, parallel pivot search and fill-in) followed by Hecke's `elementary_divisors` on the small remaining core.

Main entry point: [`sparse_elementary_divisors`](@ref). See also [`read_matrix_market`](@ref) and [`magma_elementary_divisors`](@ref) (comparison with Magma, requires a Magma installation).
"""
module SparseElementaryDivisors

using SparseArrays
import Hecke

export sparse_elementary_divisors, read_matrix_market,
       magma_elementary_divisors, timed_magma_elementary_divisors

include("solver.jl")        # core reduction
include("matrixmarket.jl")  # Matrix Market reader
include("api.jl")           # public interface
include("magma.jl")         # comparison with Magma (optional)
include("diagnostics.jl")   # consistency checks and diagnostics, not used in production

end # module
