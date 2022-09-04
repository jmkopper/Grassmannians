struct GrassmannianRing <: Generic.Ring
    k::Integer
    n::Integer
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

_valid_partition(k::Integer, n::Integer, part::Generic.Partition)::Bool = (length(part) == 0) || ((part[1] <= n-k) && (length(part) <= k))
