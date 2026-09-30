# diagnostics.jl — consistency checks and search-cost diagnostics (not part of the production path).

####################################
### Structure consistency checks ###
####################################

function check_unit_index_consistency(M::DualSparseMatrix)::Bool
    ui = M.unit_index
    ok = true
    for c ∈ 1:M.num_cols
        true_count = count(is_unit, M.col_vals[c])
        if true_count != ui.unit_count[c]
            println("Inconsistent column $c: unit_count=$(ui.unit_count[c]), actual=$true_count")
            ok = false
        end
        is_tracked = ui.buckets.bucket_of[c] != 0
        if (true_count > 0) != is_tracked
            println("Inconsistent column $c: true_count=$true_count, is_tracked=$is_tracked")
            ok = false
        end
    end
    return ok
end

function _check_buckets(db::DegreeBuckets, degrees::Vector{Int}, tracked::Vector{Bool}, name::String)::Bool
    ok = true
    for i ∈ eachindex(degrees)
        if !tracked[i]
            if db.bucket_of[i] != 0
                println("[$name] entity $i is not tracked but present in a bucket")
                ok = false
            end
            continue
        end
        bidx = db.bucket_of[i]
        if bidx != degree_to_bucket_idx(db, degrees[i])
            println("[$name] entity $i: bucket $bidx ≠ expected bucket for degree $(degrees[i])")
            ok = false
            continue
        end
        p = db.pos_in_bucket[i]
        if p < 1 || p > length(db.buckets[bidx]) || db.buckets[bidx][p] != i
            println("[$name] entity $i: inconsistent pos_in_bucket")
            ok = false
        end
    end
    return ok
end

function check_structure_consistency(M::DualSparseMatrix)::Bool
    ok = true
    nnz_rows = 0; nnz_cols = 0
    for r ∈ 1:M.num_rows
        cols = M.row_cols[r]; vals = M.row_vals[r]
        if length(cols) != length(vals)
            println("Row $r: indices/values lengths differ"); ok = false; continue
        end
        for k ∈ eachindex(cols)
            (k > 1 && cols[k - 1] >= cols[k]) && (println("Ligne $r : indices non strictement croissants"); ok = false)
            iszero(vals[k]) && (println("Row $r: stored zero value"); ok = false)
            c = cols[k]
            rr = M.col_rows[c]
            idx = searchsortedfirst(rr, r)
            if idx > length(rr) || rr[idx] != r || M.col_vals[c][idx] != vals[k]
                println("Row $r / column $c: entry not duplicated identically on the column side"); ok = false
            end
        end
        nnz_rows += length(cols)
    end
    for c ∈ 1:M.num_cols
        rr = M.col_rows[c]
        length(rr) != length(M.col_vals[c]) && (println("Column $c: indices/values lengths differ"); ok = false)
        for k ∈ 2:length(rr)
            rr[k - 1] >= rr[k] && (println("Colonne $c : indices non strictement croissants"); ok = false)
        end
        nnz_cols += length(rr)
    end
    if nnz_rows != nnz_cols
        println("row nnz ($nnz_rows) ≠ column nnz ($nnz_cols)"); ok = false
    end
    row_deg = [length(M.row_cols[r]) for r ∈ 1:M.num_rows]
    col_deg = [length(M.col_rows[c]) for c ∈ 1:M.num_cols]
    ok &= _check_buckets(M.row_buckets, row_deg, fill(true, M.num_rows), "row_buckets")
    ok &= _check_buckets(M.col_buckets, col_deg, fill(true, M.num_cols), "col_buckets")
    tracked = [M.unit_index.unit_count[c] > 0 for c ∈ 1:M.num_cols]
    ok &= _check_buckets(M.unit_index.buckets, col_deg, tracked, "unit_index")
    ok &= check_unit_index_consistency(M)
    return ok
end

"""
    state_signature(M)

Informative trajectory comparison between two instances.
"""
function state_signature(M)
    return (deepcopy(M.row_cols), deepcopy(M.row_vals), deepcopy(M.col_rows), deepcopy(M.col_vals),
            deepcopy(M.row_buckets.buckets), deepcopy(M.col_buckets.buckets),
            deepcopy(M.unit_index.buckets.buckets), copy(M.unit_index.unit_count))
end

####################################
### Pivot search cost diagnostic ###
####################################

# `diagnostic_pivot_search` reproduces EXACTLY the search (not the fill-in) of `markowitz_pivot!`; WITHOUT ANY MUTATION of `M`. `diagnose_pivot_search_cost` calls it just BEFORE the actual pivot step at each iteration: since `M` is unchanged between the two calls, both searches find the same pivot, so the mathematical result and the pivoting trajectory are EXACTLY those of a normal run with these settings; only the search is done twice (diagnostic cost ≈ 2× the search cost alone, NOT 2× the total cost). The returned NamedTuples record `n_scanned` and respect `max_buckets_scanned`.

function diagnostic_pivot_search(M::DualSparseMatrix{T};
                                  acceptable_score_threshold::Int = 4,
                                  max_buckets_scanned::Int = typemax(Int)) where {T}
    best_score::Int = typemax(Int)
    p_r::Int = -1
    p_c::Int = -1
    found_good_pivot::Bool = false
    n_scanned::Int = 0
    n_buckets_scanned::Int = 0

    max_deg::Int = M.col_buckets.max_deg
    unit_buckets::DegreeBuckets = M.unit_index.buckets

    for target_deg::Int ∈ 2:max_deg
        bidx::Int = degree_to_bucket_idx(unit_buckets, target_deg)
        bucket = unit_buckets.buckets[bidx]
        isempty(bucket) && continue
        n_buckets_scanned += 1
        for c::Int ∈ bucket
            @inbounds deg_c::Int = length(M.col_rows[c])
            deg_c != target_deg && continue
            c_rows::Vector{Int} = M.col_rows[c]
            c_vals::Vector{T} = M.col_vals[c]
            len_c_rows::Int = length(c_rows)
            for idx_r::Int ∈ 1:len_c_rows
                @inbounds val::T = c_vals[idx_r]
                if is_unit(val)
                    n_scanned += 1
                    @inbounds r::Int = c_rows[idx_r]
                    @inbounds deg_r::Int = length(M.row_cols[r])
                    score::Int = (deg_r - 1) * (deg_c - 1)
                    if score < best_score
                        best_score = score
                        p_r = r
                        p_c = c
                    end
                    if best_score <= acceptable_score_threshold
                        found_good_pivot = true
                        break
                    end
                end
            end
            found_good_pivot && break
        end
        found_good_pivot && break
        n_buckets_scanned >= max_buckets_scanned && break
    end

    early_break = found_good_pivot
    if p_r == -1
        overflow_bidx::Int = degree_to_bucket_idx(unit_buckets, max_deg + 1)
        found_good_pivot_sec::Bool = false
        for c_secours::Int ∈ unit_buckets.buckets[overflow_bidx]
            @inbounds deg_c_sec::Int = length(M.col_rows[c_secours])
            deg_c_sec <= 1 && continue
            c_rows_sec::Vector{Int} = M.col_rows[c_secours]
            c_vals_sec::Vector{T} = M.col_vals[c_secours]
            len_c_rows_sec::Int = length(c_rows_sec)
            for idx_r_sec::Int ∈ 1:len_c_rows_sec
                @inbounds val_sec::T = c_vals_sec[idx_r_sec]
                if is_unit(val_sec)
                    n_scanned += 1
                    @inbounds r_sec::Int = c_rows_sec[idx_r_sec]
                    @inbounds deg_r_sec::Int = length(M.row_cols[r_sec])
                    score_sec::Int = (deg_r_sec - 1) * (deg_c_sec - 1)
                    if score_sec < best_score
                        best_score = score_sec
                        p_r = r_sec
                        p_c = c_secours
                    end
                    if best_score <= acceptable_score_threshold
                        found_good_pivot_sec = true
                        break
                    end
                end
            end
            found_good_pivot_sec && break
        end
        early_break = found_good_pivot_sec
    end

    return (p_r = p_r, p_c = p_c, best_score = best_score,
            n_scanned = n_scanned, n_buckets_scanned = n_buckets_scanned,
            early_break = early_break, found = p_r != -1)
end

"""
    diagnose_pivot_search_cost(M; acceptable_score_threshold = 4, max_buckets_scanned = typemax(Int))

Run the complete reduction (structured elimination + Markowitz), doubling the pivot search at each iteration with `diagnostic_pivot_search` (read-only) to collect, FOR EACH MARKOWITZ PIVOT ACTUALLY APPLIED, its final score, `n_scanned`, `n_buckets_scanned`, and whether the search stopped by early break or by exhaustion (bucket or cap). MUTATES `M` — pass a copy (`deepcopy(M)`) to keep the original.

Returns a `Vector{<:NamedTuple}`, to be summarized with `summarize_pivot_search_diag`.
"""
function diagnose_pivot_search_cost(M::DualSparseMatrix{T};
                                     acceptable_score_threshold::Int = 4,
                                     max_buckets_scanned::Int = typemax(Int)) where {T}
    diag = NamedTuple{(:score, :n_scanned, :n_buckets_scanned, :early_break),
                       Tuple{Int,Int,Int,Bool}}[]
    while true
        eliminate_structured!(M)
        d = diagnostic_pivot_search(M; acceptable_score_threshold, max_buckets_scanned)
        d.found || break
        push!(diag, (score = d.best_score, n_scanned = d.n_scanned,
                     n_buckets_scanned = d.n_buckets_scanned, early_break = d.early_break))
        markowitz_pivot!(M; acceptable_score_threshold, max_buckets_scanned) || break
    end
    return diag
end

"""
    summarize_pivot_search_diag(diag; top_k=20, n_bins=8)

Summary of the results of `diagnose_pivot_search_cost`, including the statistics of `n_buckets_scanned` (in addition to `n_scanned` / `score`) — the data directly relevant to judge the effect of `max_buckets_scanned`. Prints global statistics, a Pareto-style concentration curve, a histogram by cost range, and the top-K most expensive pivots.
"""
function summarize_pivot_search_diag(diag::Vector{<:NamedTuple};
                                         top_k::Int = 20, n_bins::Int = 8)
    n_total = length(diag)
    scores = [d.score for d ∈ diag]
    scanned = [d.n_scanned for d ∈ diag]
    buckets = [d.n_buckets_scanned for d ∈ diag]
    total_scanned = sum(scanned; init = 0)
    total_buckets = sum(buckets; init = 0)

    println("Pivots applied: $n_total")
    println("Total ±1 entries scanned: $total_scanned")
    println("Total non-empty buckets scanned: $total_buckets")
    println()

    println("--- Statistiques globales ---")
    for (name, vals) ∈ (("n_scanned", scanned), ("n_buckets_scanned", buckets), ("score", scores))
        sorted_vals = sort(vals)
        n = length(sorted_vals)
        mean_v = sum(sorted_vals) / n
        median_v = sorted_vals[cld(n, 2)]
        p90 = sorted_vals[cld(90 * n, 100)]
        p99 = sorted_vals[cld(99 * n, 100)]
        println("$name: mean=$(round(mean_v, digits = 1)), median=$median_v, ",
                "p90=$p90, p99=$p99, max=$(sorted_vals[end])")
    end
    println()

    println("--- Search cost concentration (sorted by decreasing n_scanned) ---")
    order = sortperm(scanned; rev = true)
    for frac ∈ (0.001, 0.01, 0.05, 0.10, 0.25, 0.50)
        k = max(1, round(Int, frac * n_total))
        cum_k = sum(scanned[order[1:k]]; init = 0)
        println("  top $(round(100 * frac, digits = 1))% of the pivots ($k) → ",
                "$(round(100 * cum_k / max(total_scanned, 1), digits = 1))% of the total scan")
    end
    println()

    println("--- Distribution by cost range (n_scanned) ---")
    edges = [10.0^k for k ∈ 0:n_bins]
    for i ∈ 1:n_bins
        lo, hi = edges[i], edges[i + 1]
        idx = findall(x -> lo <= x < hi, scanned)
        isempty(idx) && continue
        cnt = length(idx)
        tot = sum(scanned[idx]; init = 0)
        score_range = (minimum(scores[idx]), maximum(scores[idx]))
        println("  [$(Int(lo)), $(Int(hi))[: $cnt pivots, $tot scans in total, ",
                "scores in $score_range")
    end
    println()

    println("--- The $top_k individually most expensive pivots (score, n_scanned, n_buckets_scanned) ---")
    for i ∈ 1:min(top_k, n_total)
        j = order[i]
        println("  #$i: score=$(scores[j]), n_scanned=$(scanned[j]), n_buckets_scanned=$(buckets[j])")
    end

    return (n_total = n_total, total_scanned = total_scanned, total_buckets_scanned = total_buckets)
end

"""
    sweep_search_params(M, configs)

Sweep several combinations `(acceptable_score_threshold, max_buckets_scanned)` (`configs` is a `Vector{Tuple{Int,Int}}`) on a COPY of `M` for each configuration. No `@btime`/`@benchmark` here. For each configuration, records:
  - the search cost (`total_scanned` / `total_buckets_scanned`);
  - the size of the final residual core (`core_nrows` / `core_ncols` / `core_nnz`);
  - the final divisors, to check that results agree between configurations.

Each configuration runs the diagnostic loop and the reduction inline (a single reduction per configuration) and reuses the same matrix `M_diag` to extract the core: the cost per configuration is about 2× that of the search alone.
"""
function sweep_search_params(M::DualSparseMatrix, configs::Vector{Tuple{Int,Int}})
    results = NamedTuple[]
    DiagNT = NamedTuple{(:score, :n_scanned, :n_buckets_scanned, :early_break),
                         Tuple{Int,Int,Int,Bool}}
    for (threshold, max_buckets) ∈ configs
        M_diag = deepcopy(M)
        diag = DiagNT[]
        n_ones = 0
        while true
            n_ones += eliminate_structured!(M_diag)
            d = diagnostic_pivot_search(M_diag; acceptable_score_threshold = threshold,
                                         max_buckets_scanned = max_buckets)
            d.found || break
            push!(diag, (score = d.best_score, n_scanned = d.n_scanned,
                         n_buckets_scanned = d.n_buckets_scanned, early_break = d.early_break))
            markowitz_pivot!(M_diag; acceptable_score_threshold = threshold,
                                              max_buckets_scanned = max_buckets) || break
            n_ones += 1
        end
        total_scanned = sum(d.n_scanned for d ∈ diag; init = 0)
        total_buckets = sum(d.n_buckets_scanned for d ∈ diag; init = 0)

        S_core = extract_core_smat(M_diag)
        core_nrows = Hecke.nrows(S_core)
        core_ncols = Hecke.ncols(S_core)
        core_nnz = sum(length(M_diag.row_cols[r]) for r ∈ 1:M_diag.num_rows; init = 0)
        divs = if core_nrows == 0 || core_ncols == 0
            ones(Int, n_ones)
        else
            divs_core = BigInt.(Hecke.elementary_divisors(S_core))
            divs_final = vcat(ones(Int, n_ones), divs_core)
            filter!(x -> x != 0, divs_final)
            sort!(divs_final)
        end
        push!(results, (
            threshold = threshold,
            max_buckets_scanned = max_buckets,
            n_pivots = length(diag),
            total_scanned = total_scanned,
            total_buckets_scanned = total_buckets,
            core_nrows = core_nrows,
            core_ncols = core_ncols,
            core_nnz = core_nnz,
            divisors = divs,
        ))
        println("(threshold=$threshold, max_buckets_scanned=$max_buckets) : total_scanned=$total_scanned, ",
                "total_buckets_scanned=$total_buckets, core $(core_nrows)×$(core_ncols) (nnz=$core_nnz)")
    end
    ref_divs = results[1].divisors
    for res ∈ results
        if res.divisors != ref_divs
            println("WARNING: results differ for (threshold=$(res.threshold), ",
                    "max_buckets_scanned=$(res.max_buckets_scanned)) — expected $ref_divs, ",
                    "got $(res.divisors)")
        end
    end
    return results
end

# NOT INCLUDED (to be added if really needed): the fill-in size collection and summary tools. To judge the fill-in quality of a configuration (threshold, max_buckets_scanned), the size of the final residual core (core_nrows / core_ncols / core_nnz, already reported by `sweep_search_params`) is enough.

