# magma.jl — comparison with Magma (optional: requires a Magma installation, local or through SSH).

# Magma's "flat" sparse format: for each row, the number of non-zero entries followed by the (column, value) pairs.
function _magma_matrix_definition(rows::RowsData, m::Int, n::Int)
    Q = Int[]
    for row ∈ rows
        push!(Q, length(row))
        for (c, v) ∈ row
            push!(Q, c); push!(Q, v)
        end
    end
    return "M := SparseMatrix($(m), $(n), [ " * join(Q, ", ") * " ]);"
end

function _parse_magma_output(raw::AbstractString)
    m_err = match(r"##MAGMA_ERROR_BEGIN##(.*?)##MAGMA_ERROR_END##"s, raw)
    m_err === nothing || error("Magma reported an error during the computation:\n" * strip(m_err.captures[1]))
    t_exec = NaN
    m_time = match(r"Time:\s*([0-9.eE+\-]+)", raw)
    m_time === nothing || (t_exec = parse(Float64, m_time.captures[1]))
    divs = BigInt[]
    m_divs = match(r"##DIVS_BEGIN##(.*?)##DIVS_END##"s, raw)
    if m_divs !== nothing
        cleaned = strip(replace(m_divs.captures[1], r"[\[\]\r\n]" => ""))
        if !isempty(cleaned)
            occursin(r"^[\s,0-9+\-]+$", cleaned) ||
                error("non-numeric output from Magma:\n" * raw)
            divs = parse.(BigInt, strip.(split(cleaned, ",")))
        end
    end
    return divs, t_exec
end

"""
    timed_magma_elementary_divisors(A; ssh_host = nothing, magma_command = "magma -b",
                                    connect_timeout = 15)
    timed_magma_elementary_divisors(rows, m, n; kwargs...)

Compute the elementary divisors of `A` with Magma, to compare results and timings. Return `(divisors::Vector{BigInt}, t_magma::Float64)`, where `t_magma` is the time measured BY MAGMA for the single call to `ElementaryDivisors` (excluding the connection, the start-up of Magma and the transfer of the matrix; `NaN` if it could not be read). `A` can be of any type accepted by [`sparse_elementary_divisors`](@ref), including a Matrix Market file name.

Without `ssh_host`, `magma_command` is run locally; otherwise the command is run through `ssh -o ConnectTimeout=… ssh_host magma_command`. The matrix is sent on the standard input.

See [`magma_elementary_divisors`](@ref) for a version returning only the divisors.
"""
function timed_magma_elementary_divisors(rows::RowsData, m::Int, n::Int;
                                         ssh_host::Union{Nothing,AbstractString} = nothing,
                                         magma_command::AbstractString = "magma -b",
                                         connect_timeout::Int = 15)
    script = """
    $(_magma_matrix_definition(rows, m, n))
    try
        time divs := ElementaryDivisors(M);
        print "##DIVS_BEGIN##";
        print divs;
        print "##DIVS_END##";
    catch e
        print "##MAGMA_ERROR_BEGIN##";
        print e;
        print "##MAGMA_ERROR_END##";
    end try;
    quit;
    """
    cmd = ssh_host === nothing ? Cmd(String.(split(magma_command))) :
          Cmd(["ssh", "-o", "ConnectTimeout=$(connect_timeout)", String(ssh_host), String(magma_command)])
    out, err = IOBuffer(), IOBuffer()
    proc = nothing
    process_error = nothing
    try
        # `ignorestatus`: Magma may return a non-zero exit code even though the error was caught
        proc = run(pipeline(ignorestatus(cmd); stdin = IOBuffer(script), stdout = out, stderr = err))
    catch e
        process_error = e
    end
    result = String(take!(out))
    errmsg = String(take!(err))
    usable = occursin("##MAGMA_ERROR_BEGIN##", result) || occursin("##DIVS_BEGIN##", result)
    if process_error !== nothing || (!usable && proc !== nothing && !success(proc))
        error("""
        The call to Magma failed (command: $cmd)
        --- stderr ---
        $(isempty(strip(errmsg)) ? "[empty]" : errmsg)
        --- stdout ---
        $(isempty(strip(result)) ? "[empty]" : result)
        Exception: $process_error
        (with SSH, test: ssh $(ssh_host === nothing ? "<host>" : ssh_host) "echo ok")
        """)
    end
    return _parse_magma_output(result)
end

timed_magma_elementary_divisors(m::Int, n::Int, rows::RowsData; kwargs...) =
    timed_magma_elementary_divisors(rows, m, n; kwargs...)
timed_magma_elementary_divisors(A::SparseMatrixCSC{<:Integer}; kwargs...) =
    timed_magma_elementary_divisors(_rows_data(A), size(A)...; kwargs...)
timed_magma_elementary_divisors(A::AbstractMatrix{<:Integer}; kwargs...) =
    timed_magma_elementary_divisors(sparse(A); kwargs...)
timed_magma_elementary_divisors(S::Union{Hecke.SMat{Hecke.ZZRingElem},Hecke.ZZMatrix}; kwargs...) =
    timed_magma_elementary_divisors(_rows_data(S), Hecke.nrows(S), Hecke.ncols(S); kwargs...)
timed_magma_elementary_divisors(path::AbstractString; kwargs...) =
    timed_magma_elementary_divisors(read_matrix_market(path); kwargs...)

"""
    magma_elementary_divisors(A; verbose = false, kwargs...)

Non-zero elementary divisors of `A` computed by Magma, as a `Vector{BigInt}` (same result type as [`sparse_elementary_divisors`](@ref)). With `verbose = true`, Magma's own computation time is printed. The other keyword arguments (`ssh_host`, `magma_command`, `connect_timeout`) and the accepted inputs are those of [`timed_magma_elementary_divisors`](@ref), which also returns the time.
"""
function magma_elementary_divisors(args...; verbose::Bool = false, kwargs...)
    divs, t = timed_magma_elementary_divisors(args...; kwargs...)
    verbose && println("Magma: ElementaryDivisors computed in ", t, " s")
    return divs
end
