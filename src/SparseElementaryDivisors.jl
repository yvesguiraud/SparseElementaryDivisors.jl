"""
    SparseElementaryDivisors

Elementary divisors (Smith normal form invariants) of large **sparse integer matrices**, by exact sparse elimination (unit pivots with Markowitz-style ordering, parallel pivot search and fill-in) followed by Hecke's `elementary_divisors` on the small remaining core.

References (details in the README): Markowitz (1957) for the pivot score; LaMacchia and Odlyzko (1991) for the structured elimination; Dumas, Saunders and Villard (2001) on the Smith normal form of large sparse integer matrices; Fieker, Hart, Hofmann and Johansson (2017) for Hecke, with FLINT.

Main entry point: [`sparse_elementary_divisors`](@ref), with its settings [`ReductionConfig`](@ref) and its statistics [`ReductionStats`](@ref). See also [`read_matrix_market`](@ref), [`write_matrix_market`](@ref) and [`magma_elementary_divisors`](@ref) (comparison with Magma, requires a Magma installation).
"""
module SparseElementaryDivisors

using SparseArrays
import Hecke

export sparse_elementary_divisors, ReductionConfig, ReductionStats, CoefficientOverflowError,
       read_matrix_market, write_matrix_market, magma_elementary_divisors,
       timed_magma_elementary_divisors

# Public without being exported (`public` exists from Julia 1.11 on; the package supports 1.10,
# where the word would be a syntax error, hence the expression built and evaluated here).
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Expr(:public, :DEFAULT_VALUE_TYPES, :ReductionTrace, :RowsData))
end

include("errors.jl")        # dedicated error (overflow of the coefficients)
include("solver.jl")        # core reduction
include("matrixmarket.jl")  # Matrix Market reader and writer
include("api.jl")           # public interface
include("magma.jl")         # comparison with Magma (optional)
include("diagnostics.jl")   # consistency checks and diagnostics, not used in production

end # module
