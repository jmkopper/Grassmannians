#========================
    Vector operations
========================#
function ==(a::SchubertCycle, b::SchubertCycle)
    parent(a) != parent(b) && return false
        for (p, v) in a.terms
            if haskey(b.terms, p)
                b.terms[p] == v || return false
            elseif v !=0
                return false
            end
        end
    return true
end


function +(a::SchubertCycle, b::SchubertCycle)::SchubertCycle
    c = parent(a)()
    for (part, v) in a.terms
        if haskey(b.terms, part)
            c.terms[part] = v + b.terms[part]
        else
            c.terms[part] = v
        end
    end
    for (part,v) in b.terms
        if !haskey(c.terms, part)
            c.terms[part] = v
        end
    end
    return c
end

function *(a::Integer, b::SchubertCycle)::SchubertCycle
    terms = Dict(p=>a*v for (p,v) in b.terms)
    return SchubertCycle(b.k, b.n, terms, parent(b))
end
*(a::SchubertCycle, b::Integer)::SchubertCycle = *(b,a)
-(a::SchubertCycle) = (-1)*a
-(a::SchubertCycle, b::SchubertCycle)::SchubertCycle = a + (-1)*b

#========================
    Schubert calculus
========================#

#= Pieri's formula =#
@inline _is_pieri(p::Generic.Partition)::Bool = (length(p) <= 1)

function _valid_pieri_summand(product_partition::Generic.Partition, summand_partition::Generic.Partition)::Bool
    """
    Determine whether `summand_partition` is a valid summand in Pieri's formula, assuming its codimension is correct (not checked)
    """
    for i in 1:max(length(product_partition), length(summand_partition))
        cur_prod = length(product_partition) >= i ? product_partition[i] : 0
        prev_prod = (length(product_partition) >= i-1) && (i > 1) ? product_partition[i-1] : 0
        cur_sum = length(summand_partition) >= i ? summand_partition[i] : 0
        if cur_sum < cur_prod || (cur_sum > prev_prod && i > 1)
            return false
        end
    end
    return true
end

function _pieri_prod(p::Generic.Partition, q::Generic.Partition)::Vector{Generic.Partition}
    v_n = p.n + q.n
    valid_partitions = []
    if _is_pieri(q)
        p, q = q, p
    end
    for part in Generic.partitions(v_n)
        if _valid_pieri_summand(q, part)
            push!(valid_partitions, part)
        end
    end

    return valid_partitions
end

#= Giambelli's formula =#

function _giambelli_matrix(g::GrassmannianRing, p::Generic.Partition)
    mat_list = []
    k = length(p)
    for i in 1:k
        row = []
        for j in 1:k
            l = p[i]
            if l-i+j > 0 && l-i+j <= g.n-g.k
                push!(row, g([l-i+j]))
            elseif l-i+j == 0
                push!(row, one(g))
            else
                push!(row, zero(g))
            end
        end
        push!(mat_list, row)
    end
    return matrix(g, reduce(vcat, transpose.(mat_list)))
end

function _get_cofactor(R::Generic.Ring, M, row::Integer, col::Integer)
    n = size(M)[1]
    cofactor = []

    for i in 1:n
        if i != row
            new_row = []
            for j in 1:n
                if j != col
                    push!(new_row, M[i,j])
                end
            end
            push!(cofactor, new_row)
        end
    end

    return matrix(R, reduce(vcat, transpose.(cofactor)))
end


function _cofactor_det(R::Generic.Ring, M)
    n = size(M)[1]
    if n == 1
        return M[1,1]
    end

    det = 0
    for col in 1:n
        N = _get_cofactor(R, M, 1, col)
        N[1, :] = N[1, :] * M[1, col]
        det += _cofactor_det(R, N) * (-1)^(1+col)
    end
    return det
end

@inline function _giambellify(g::GrassmannianRing, p::Generic.Partition, q::Generic.Partition)
    m = _giambelli_matrix(g, p)
    m[1, :] = m[1, :]*g(q)
    return m
end

@inline function *(p::Generic.Partition, q::Generic.Partition)::Vector{Generic.Partition}
    if !_is_pieri(p) && !_is_pieri(q)
        print("requires a special partition")
        return
    end

    return _pieri_prod(p, q)
end

#=
    Product formula for Schubert cycles. Implements Giambelli and Pieri formulas.
    Cofactor expansion algorithm for determinants is used to force Pieri's rule in
    the actual computation of Giambelli's formula.
=#
function *(a::SchubertCycle, b::SchubertCycle)::SchubertCycle
    if parent(a) != parent(b)
        throw(DomainError(a, "cycles must be on the same Grassmannian"))
    end

    g = parent(a)

    if a == zero(g) || b == zero(g)
        return zero(g)
    end

    terms = Dict{Generic.Partition, Integer}()
    c = g(0)

    for (at, bt) in Iterators.product(_nonzero_terms(a), _nonzero_terms(b))
        (part_a, coeff_a) = at
        (part_b, coeff_b) = bt

        if _is_pieri(part_a) || _is_pieri(part_b)
            for term in _pieri_prod(part_a, part_b)
                if _valid_partition(g.k, g.n, term)
                    terms[term] = haskey(terms, term) ? terms[term] + coeff_a*coeff_b : coeff_a*coeff_b
                end
            end

        else
            c += _cofactor_det(g, _giambellify(g, part_a, part_b))
        end

    end

    d = SchubertCycle(a.k, a.n, terms, parent(a))
    x = c + d
    _remove_zero_terms!(x)
    return x
end

function ^(a::SchubertCycle, n::Integer)::SchubertCycle
    n < 0 && throw(DomainError(n, "exponent must be positive integer"))
    prod = one(parent(a))
    for _ in 1:n
        prod *= a
    end
    return prod
end


