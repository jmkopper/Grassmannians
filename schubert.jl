import AbstractAlgebra.parent
import AbstractAlgebra: zero, one, elem_type, parent_type

struct GrassmannianRing <: Generic.Ring
    k::Integer
    n::Integer

    function GrassmannianRing(k, n)
        # valid_cycles = Dict(p => 0 for p in all_valid_partitions(k, n))
        new(k, n)
    end
end


struct SchubertCycle <: RingElem
    k::Integer
    n::Integer
    terms::Dict{Generic.Partition, Integer}
    parent::GrassmannianRing
end

###############################################################################
#
#   Implementing Ring/RingElem interface for Grassmannians
#
###############################################################################

elem_type(::Type{GrassmannianRing}) = SchubertCycle
parent_type(::Type{SchubertCycle}) = GrassmannianRing

parent(a::SchubertCycle) = a.parent

(g::GrassmannianRing)(x::Integer)::SchubertCycle = SchubertCycle(g.k, g.n, Dict(Generic.Partition([0], false)=>x), g)
one(g::GrassmannianRing)::SchubertCycle = g(1)
zero(g::GrassmannianRing)::SchubertCycle = g(0)
is_zero(a::SchubertCycle) = all(values(a.terms) .== 0)
(g::GrassmannianRing)() = zero(g)


function (g::GrassmannianRing)(p::Generic.Partition)::SchubertCycle
    if length(p) <= g.k && p[1] <= g.n-g.k
        return SchubertCycle(g.k, g.n, Dict(p => 1), g)
    else
        throw(DomainError(g, "Invalid partition $p for $g"))
    end
end

(g::GrassmannianRing)(p::Vector)::SchubertCycle = g(Partition(p))


function (g::GrassmannianRing)(a::SchubertCycle)::SchubertCycle
    parent(a) == g && return a
    throw(DomainError(g, "Schubert cycle $a is not an element of $g"))
end

#=
Thots:
- Have the grassmannian constructor initialize ALL possible partitions as a hash with k=>0 for all k
- Change *(...) to iterate through this instead of re-generate all partitions
- Have schubertcycle constructor pull from this too
- use Holy traits to detect special partitions for pieri
=#


###############################################################################
#
#   Vector operations
#
###############################################################################
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

###############################################################################
#
#   Schubert calculus
#
###############################################################################

function is_special_cycle(a::SchubertCycle)::Bool
    nonzero_terms = [v for v in values(a.terms) if v != 0]
    return length(nonzero_terms) == 1
end

_is_pieri(p::Generic.Partition)::Bool = (length(p) <= 1)

   
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

function _giambelli_matrix(g::GrassmannianRing, p::Generic.Partition)
    mat_list = []
    k = length(p)
    for i in 1:k
        row = []
        for j in 1:k
            l = p[i]
            if l-i+j > 0 && l-i+j <= g.n-g.k
                push!(row, g([l-i+j]))
            else
                push!(row, one(g))
            end
        end
        push!(mat_list, row)
    end
    return matrix(g, reduce(vcat, transpose.(mat_list)))
end

function _giambellify(g::GrassmannianRing, p::Generic.Partition, q::Generic.Partition)
    m = _giambelli_matrix(g, p)
    m[1, :] = m[1, :]*g(q)
    return m
end

function *(p::Generic.Partition, q::Generic.Partition)::Vector{Generic.Partition}
    if !_is_pieri(p) && !_is_pieri(q)
        print("Requires a special partition")
        return
    end
    if !_is_pieri(p)
        # Make p the special partition
        p, q = q, p
    end

    return _pieri_prod(p, q)
end

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
    for (part_a, coeff_a) in a.terms
        for (part_b, coeff_b) in b.terms
            if _is_pieri(part_a) || _is_pieri(part_b)
                prod = _pieri_prod(part_a, part_b)
                for term in prod
                    if _valid_partition(g.k, g.n, term)
                        terms[term] = haskey(terms, term) ? terms[term] + coeff_a*coeff_b : coeff_a*coeff_b
                    end
                end
            else
                c += det(_giambellify(g, part_a, part_b))
            end
        end
    end

    d = SchubertCycle(a.k, a.n, terms, parent(a))
    return c+ d
end

# Pretty printing of SchubertCycle
function Base.show(io::IO, a::SchubertCycle)
    if is_zero(a)
        print(io, "0")
        return
    end
    str::String = ""
    for (term, coeff) in a.terms
        term_str::String = "σ("
        for (index, val) in enumerate(term)
            term_str *= "$val"
            if index < length(term)
                term_str *= ","
            end
        end
        term_str *= ")"
        if coeff == 1
            str *= term_str * " + "
        elseif coeff == -1
            str *= "-" * term_str * " + "
        elseif coeff != 0
            str *= string(coeff) * term_str * " + "
        end
    end
    print(io, str[1:end-2])
end

function test(a::SchubertCycle)::Integer
    print(a.partition)
    return a.k
end

_valid_partition(k::Integer, n::Integer, part::Generic.Partition)::Bool = (length(part) == 0) || ((part[1] <= n-k) && (length(part) <= k))
_valid_schubert_cylce(k::Integer, n::Integer, a::SchubertCycle)::Bool = all(_valid_partition(k, n, part) for part in a.terms)

function all_valid_partitions(k::Integer, n::Integer)
    parts = []
    for i in 1:k*(n-k)
        for p in Generic.partitions(i)
            _valid_partition(k, n, p) && push!(parts, p)
        end
    end
    push!(parts, Partition([0], false))
    return parts
end