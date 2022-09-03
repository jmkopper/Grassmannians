import AbstractAlgebra.parent

struct GrassmannianRing <: Generic.Ring
    k::Integer
    n::Integer

    GrassmannianRing(k, n) = new(k, n)
end


struct SchubertCycle <: RingElem
    k::Integer
    n::Integer
    terms::Dict{Generic.Partition, Integer}
    parent::GrassmannianRing

    # SchubertCycle(k, n, terms, parent) = new(k, n, Dict(sort(collect(terms), by = x->x[1])), parent)
end

###############################################################################
#
#   Implementing Ring/RingElem interface for Grassmannians
#
###############################################################################

elem_type(::Type{GrassmannianRing}) = SchubertCycle
parent_type(::Type{SchubertCycle}) = GrassmannianRing

parent(a::SchubertCycle) = a.parent

function zero(g::GrassmannianRing)
    qq = Dict{Generic.Partition, Integer}()
    return SchubertCycle(g.k, g.n, qq, g)
end

function one(g::GrassmannianRing)
    p = Partition([0], false)
    return SchubertCycle(g.k, g.n, Dict(p=> 1), g)
end

function is_zero(a::SchubertCycle)
    for (_, v) in a.terms
        if v > 0
            return false
        end
    end
    return true
end

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
    if parent(a) == g
        return a
    end
    throw(DomainError(g, "Schubert cycle $a is not an element of $g"))
end

function ==(a::SchubertCycle, b::SchubertCycle)
    if parent(a) != parent(b)
        return false
    elseif a.terms != b.terms
        return false
    else
        return true
    end
end

###############################################################################
#
#   Schubert calculus
#
###############################################################################

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


function *(a::SchubertCycle, b::SchubertCycle)::SchubertCycle
    if parent(a) != parent(b)
        throw(DomainError(a, "cycles must be on the same Grassmannian"))
    end

    g = parent(a)

    if a == zero(g)
        return b
    end
    if b == zero(g)
        return a
    end

    terms = Dict{Generic.Partition, Integer}()
    for (part_a, coeff_a) in a.terms
        for (part_b, coeff_b) in b.terms
            for term in part_a*part_b
                if _valid_partition(a.k, a.n, term)
                    terms[term] = haskey(terms, term) ? terms[term] + coeff_a*coeff_b : coeff_a*coeff_b
                end
            end
        end
    end
    return SchubertCycle(a.k, a.n, terms, parent(a))
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
        str *= coeff != 1 ? string(coeff) * term_str * " + " : term_str * " + "
    end
    print(io, str[1:end-2])
end

function test(a::SchubertCycle)::Integer
    print(a.partition)
    return a.k
end

_valid_partition(k::Integer, n::Integer, part::Generic.Partition)::Bool = (part[1] <= n-k) && (length(part) <= k)
_valid_schubert_cylce(k::Integer, n::Integer, a::SchubertCycle)::Bool = all(_valid_partition(part) for part in terms)