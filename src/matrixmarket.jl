# matrixmarket.jl — reading integer matrices in Matrix Market format (coordinate or array).

"""
    read_matrix_market(path_or_io) -> SparseMatrixCSC{Int,Int}

Read a matrix in [Matrix Market](https://math.nist.gov/MatrixMarket/formats.html) format: `coordinate` (the format suited to sparse matrices) or `array`, of type `integer` or `pattern` (all entries equal to 1); `real` is accepted if all values are integers. The symmetry can be `general`, `symmetric` or `skew-symmetric` (the latter two for `coordinate` only). Comment lines (`%`) are ignored.
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
