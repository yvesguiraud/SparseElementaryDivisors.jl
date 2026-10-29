# matrixmarket.jl — reading integer matrices in Matrix Market format (coordinate or array),
# and writing them (coordinate integer general).

"""
    read_matrix_market(path_or_io) -> SparseMatrixCSC{Int,Int}

Read a matrix in [Matrix Market](https://math.nist.gov/MatrixMarket/formats.html) format (Boisvert, Pozo and Remington 1996, references in the README): `coordinate` (the format suited to sparse matrices) or `array`, of type `integer` or `pattern` (all entries equal to 1); `real` is accepted if all values are integers. The symmetry can be `general`, `symmetric` or `skew-symmetric` (the latter two for `coordinate` only). Comment lines (`%`) are ignored. The values are read as `Int`: an entry beyond a machine word raises an error. See also [`write_matrix_market`](@ref).
"""
function read_matrix_market(io::IO)
    header = readline(io)
    startswith(header, "%%MatrixMarket") ||
        throw(ArgumentError("missing Matrix Market header (first line: \"$header\")"))
    tok = lowercase.(split(header))
    length(tok) >= 5 || throw(ArgumentError("incomplete Matrix Market header: \"$header\""))
    obj, fmt, field, sym = tok[2], tok[3], tok[4], tok[5]
    obj == "matrix" || throw(ArgumentError("unsupported object \"$obj\" (only \"matrix\")"))
    field ∈ ("integer", "real", "pattern") ||
        throw(ArgumentError("unsupported field \"$field\" (integer, pattern, real)"))
    sym ∈ ("general", "symmetric", "skew-symmetric") ||
        throw(ArgumentError("unsupported symmetry \"$sym\""))

    _value(s) = begin
        field == "pattern" && return 1
        field == "integer" && return parse(Int, s)
        x = parse(Float64, s)
        isinteger(x) || throw(ArgumentError("non-integer value \"$s\" in a \"real\" file"))
        return Int(x)
    end

    # first non-comment, non-empty line: dimensions
    size_line = ""
    for line in eachline(io)
        (isempty(strip(line)) || startswith(line, "%")) && continue
        size_line = line
        break
    end
    isempty(size_line) && throw(ArgumentError("missing dimensions line"))
    dims = parse.(Int, split(size_line))

    if fmt == "coordinate"
        length(dims) == 3 || throw(ArgumentError("invalid dimensions line: \"$size_line\""))
        m, n, nz = dims
        I = Int[]; J = Int[]; V = Int[]
        sizehint!(I, sym == "general" ? nz : 2nz); sizehint!(J, sym == "general" ? nz : 2nz)
        sizehint!(V, sym == "general" ? nz : 2nz)
        count = 0
        for line in eachline(io)
            (isempty(strip(line)) || startswith(line, "%")) && continue
            f = split(line)
            i = parse(Int, f[1]); j = parse(Int, f[2])
            v = field == "pattern" ? 1 : _value(f[3])
            (1 <= i <= m && 1 <= j <= n) || throw(ArgumentError("index out of range: \"$line\""))
            push!(I, i); push!(J, j); push!(V, v)
            if sym != "general" && i != j
                push!(I, j); push!(J, i); push!(V, sym == "symmetric" ? v : -v)
            end
            count += 1
        end
        count == nz || throw(ArgumentError("$count entries read, $nz announced"))
        return sparse(I, J, V, m, n)
    elseif fmt == "array"
        sym == "general" || throw(ArgumentError("array format: only the \"general\" symmetry is supported"))
        length(dims) == 2 || throw(ArgumentError("invalid dimensions line: \"$size_line\""))
        m, n = dims
        vals = Int[]
        sizehint!(vals, m * n)
        for line in eachline(io)
            (isempty(strip(line)) || startswith(line, "%")) && continue
            push!(vals, _value(strip(line)))
        end
        length(vals) == m * n || throw(ArgumentError("$(length(vals)) values read, $(m * n) expected"))
        return sparse(reshape(vals, m, n))          # column-major order
    else
        throw(ArgumentError("unsupported format \"$fmt\" (coordinate, array)"))
    end
end

read_matrix_market(path::AbstractString) = open(read_matrix_market, path)

"""
    write_matrix_market(path::AbstractString, A; comment = nothing)
    write_matrix_market(io::IO, A; comment = nothing)

Write the integer matrix `A` (a `SparseMatrixCSC{<:Integer}`, or any `AbstractMatrix{<:Integer}`,
converted by `sparse`) in [Matrix Market](https://math.nist.gov/MatrixMarket/formats.html)
format (Boisvert, Pozo and Remington 1996), `coordinate integer general`: the header line, the comment lines, the line `m n nnz`,
then one line `i j v` per nonzero entry, sorted by row, then by column. The stored zeros are
left out (`nnz` counts the entries written); the values are written in decimal, whatever their
integer type (`BigInt` included); the lines end with `\\n` on every system. The same matrix
written with the same comment always gives the same bytes.

`comment`: `nothing`, a string (several lines if it contains `\\n`) or a vector of strings; each
line is written after `% `. An existing file is overwritten. Returns `nothing`.

[`read_matrix_market`](@ref) reads such a file back, with values of type `Int`: an entry beyond
a machine word can be written but not read back yet.

```julia
write_matrix_market("d2.mtx", A; comment = "d_2 of B*(G34), trivial coefficients")
```
"""
function write_matrix_market(io::IO, A::SparseMatrixCSC{<:Integer}; comment = nothing)
    lines = _comment_lines(comment)
    m, n = size(A)
    At = copy(transpose(A))        # its columns are the rows of A, column indices sorted
    rv = rowvals(At)
    nz = nonzeros(At)
    print(io, "%%MatrixMarket matrix coordinate integer general\n")
    for l ∈ lines
        print(io, "% ", l, "\n")
    end
    print(io, m, " ", n, " ", count(!iszero, nz), "\n")
    for i ∈ 1:m, k ∈ nzrange(At, i)
        v = nz[k]
        iszero(v) && continue
        print(io, i, " ", rv[k], " ", _mm_value(v), "\n")
    end
    return nothing
end

write_matrix_market(io::IO, A::AbstractMatrix{<:Integer}; comment = nothing) =
    write_matrix_market(io, sparse(A); comment)

function write_matrix_market(path::AbstractString, A; comment = nothing)
    open(path, "w") do io
        write_matrix_market(io, A; comment)
    end
    return nothing
end

# A Bool is written as 0 or 1, not as false or true.
_mm_value(v::Integer) = v
_mm_value(v::Bool) = Int(v)

# The lines of the comment: a string is split at its line breaks (a final one is dropped).
function _comment_lines(comment)
    lines = String[]
    comment === nothing && return lines
    for part ∈ (comment isa AbstractString ? (comment,) : comment)
        part isa AbstractString ||
            throw(ArgumentError("comment must be nothing, a string or a vector of strings (got $(repr(comment)))"))
        for l ∈ split(chomp(part), '\n')
            push!(lines, rstrip(l, '\r'))
        end
    end
    return lines
end
