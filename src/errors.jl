# errors.jl — the dedicated error of the solver.

"""
    CoefficientOverflowError(value_type, unit_pivots, entry)

The integer coefficients do not fit in `value_type`. Raised by
[`sparse_elementary_divisors`](@ref) when the last type of the sequence given by its keyword
`value_type` does not suffice; with the default `value_type = :auto`, whose sequence ends
with `BigInt`, it is never raised.

- `value_type`: the type that overflowed;
- `unit_pivots`: number of unit pivots done before the overflow (`0` for an input entry);
- `entry`: `nothing` for an overflow during the elimination, or `(row, column, value)` of
  an input entry that does not fit in `value_type`.
"""
struct CoefficientOverflowError <: Exception
    value_type::Type
    unit_pivots::Int
    entry::Union{Nothing, Tuple{Int, Int, BigInt}}
end

# The message starts as that of v0.1.0 ("Arithmetic overflow detected with value_type=…"),
# which scripts may recognise.
function Base.showerror(io::IO, e::CoefficientOverflowError)
    print(io, "CoefficientOverflowError: Arithmetic overflow detected with value_type=", e.value_type)
    if e.entry === nothing
        print(io, ", after ", e.unit_pivots, " unit pivots.")
    else
        r, c, v = e.entry
        print(io, ": input entry (", r, ", ", c, ") = ", v, " does not fit.")
    end
    print(io, " Retry with a wider value_type, or with value_type = :auto (the default), ",
          "which widens automatically.")
end
