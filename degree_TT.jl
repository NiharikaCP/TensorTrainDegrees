
#----------------------
# k-subset generator 
""""
All size-`k` subsets of `1:n`, each returned as a sorted Vector{Int},
in lexicographic order. 
"""

function k_subsets(n::Int, k::Int)
    result = Vector{Vector{Int}}()
    if k < 0 || k > n
        return result
    end
    comb = collect(1:k)
    while true
        push!(result, copy(comb))
        i = k
        while i >= 1 && comb[i] == n - k + i
            i -= 1
        end
        i < 1 && break
        comb[i] += 1
        for j in (i + 1):k
            comb[j] = comb[j - 1] + 1
        end
    end
    return result
end

#----------------------
# Lattice paths

"""
All lattice paths from cell (r0,c0) to cell (r1,c1) using unit "up"
(row -= 1) or "right" (col += 1) steps. Requires r0 ≥ r1 and c1 ≥ c0.
Returns a Vector of paths, each Vector{Tuple{Int,Int}} of the
cells visited.
"""
function lattice_paths(r0::Int, c0::Int, r1::Int, c1::Int)
    up_steps    = r0 - r1
    right_steps = c1 - c0
    @assert up_steps >= 0 && right_steps >= 0 "invalid path endpoints"
    n = up_steps + right_steps
    paths = Vector{Vector{Tuple{Int,Int}}}()
    for rightpos in k_subsets(n, right_steps)
        rightset = Set(rightpos)
        cells = Vector{Tuple{Int,Int}}(undef, n + 1)
        r, c = r0, c0
        cells[1] = (r, c)
        for step in 1:n
            if step in rightset
                c += 1
            else
                r -= 1
            end
            cells[step + 1] = (r, c)
        end
        push!(paths, cells)
    end
    return paths
end

#-------------------------------
# Non-intersecting path families

"""
All t-tuples of pairwise non-intersecting (no shared cell) lattice paths
(P_1,...,P_t) on a p×q matrix, where P_ℓ connects the ℓ-th entry of the
last row, (p,ℓ), to the ℓ-th entry of the last column, (ℓ,q), for
ℓ = 1,...,t. Requires t ≤ min(p,q). 

Returns a Vector of families, each family a Vector{Vector{Tuple{Int,Int}}}
of the t paths.
"""
function noncrossing_path_families(p::Int, q::Int, t::Int)
    @assert t <= min(p, q) "need t ≤ min(p,q)"
    candidate_paths = [lattice_paths(p, ell, ell, q) for ell in 1:t]

    results = Vector{Vector{Vector{Tuple{Int,Int}}}}()
    current = Vector{Vector{Tuple{Int,Int}}}()
    used = Set{Tuple{Int,Int}}()

    function backtrack(ell::Int)
        if ell > t
            push!(results, deepcopy(current))
            return
        end
        for path in candidate_paths[ell]
            pathset = Set(path)
            if isempty(intersect(pathset, used))
                push!(current, path)
                union!(used, pathset)
                backtrack(ell + 1)
                setdiff!(used, pathset)
                pop!(current)
            end
        end
    end

    backtrack(1)
    return results
end


# Mapping matrix cells back to index triples (j1,j2,j3)

function cell_to_tensor_psi1(cell::Tuple{Int,Int}, d2::Int, d3::Int)
    j1, c = cell
    j2 = div(c - 1, d3) + 1
    j3 = mod(c - 1, d3) + 1
    return (j1, j2, j3)
end

function cell_to_tensor_psi2(cell::Tuple{Int,Int}, d1::Int)
    r, j3 = cell
    j2 = div(r - 1, d1) + 1
    j1 = mod(r - 1, d1) + 1
    return (j1, j2, j3)
end

#--------------------
# Degree computation

"""
Affine dimension of the order-3 tensor train variety TT_{(d1,d2,d3),(r1,r2)} 
"""
tt_dimension(d1, d2, d3, r1, r2) = d1*r1 + r1*r2*d2 + r2*d3 - r1^2 - r2^2

"""
Degree of the order-3 tensor train variety TT_{(d1,d2,d3),(r1,r2)}, via
the combinatorial method presented in this paper
"""
function compute_degree(d1::Int, d2::Int, d3::Int, r1::Int, r2::Int; verbose::Bool=false)
    p1, q1 = d1, d2 * d3
    p2, q2 = d1 * d2, d3
    @assert r1 <= min(p1, q1) "r1 out of range"
    @assert r2 <= min(p2, q2) "r2 out of range"

    fam1 = noncrossing_path_families(p1, q1, r1)
    fam2 = noncrossing_path_families(p2, q2, r2)
    dim  = tt_dimension(d1, d2, d3, r1, r2)

    if verbose
        println("#families on ψ^(1): ", length(fam1))
        println("#families on ψ^(2): ", length(fam2))
        println("target dimension:   ", dim)
    end

    # Precompute, for each family, the set of tensor-index triples it covers.
    fam1_sets = [Set(cell_to_tensor_psi1(cell, d2, d3) for path in fam for cell in path)
                 for fam in fam1]
    fam2_sets = [Set(cell_to_tensor_psi2(cell, d1) for path in fam for cell in path)
                 for fam in fam2]

    deg = 0
    for s1 in fam1_sets, s2 in fam2_sets
        if length(intersect(s1, s2)) == dim
            deg += 1
        end
    end
    return deg
end

