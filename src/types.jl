#=
    Implementing Ring/RingElem interface for Grassmannians
=#

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


elem_type(::Type{GrassmannianRing}) = SchubertCycle
parent_type(::Type{SchubertCycle}) = GrassmannianRing
parent(a::SchubertCycle) = a.parent

(g::GrassmannianRing)(x::Integer)::SchubertCycle = SchubertCycle(g.k, g.n, Dict(Generic.Partition([0], false)=>x), g)
(g::GrassmannianRing)() = zero(g)

one(g::GrassmannianRing)::SchubertCycle = g(1)
zero(g::GrassmannianRing)::SchubertCycle = g(0)
is_zero(a::SchubertCycle) = all(values(a.terms) .== 0)

(g::GrassmannianRing)(p::Vector)::SchubertCycle = g(Partition(p))
function (g::GrassmannianRing)(p::Generic.Partition)::SchubertCycle
    if length(p) <= g.k && p[1] <= g.n-g.k
        return SchubertCycle(g.k, g.n, Dict(p => 1), g)
    else
        throw(DomainError(g, "Invalid partition $p for $g"))
    end
end


function (g::GrassmannianRing)(a::SchubertCycle)::SchubertCycle
    parent(a) == g && return a
    throw(DomainError(g, "Schubert cycle $a is not an element of $g"))
end

@inline _valid_partition(k::Integer, n::Integer, part::Generic.Partition)::Bool = (length(part) == 0) || ((part[1] <= n-k) && (length(part) <= k))

#=
    Geometry
=#


@inline dim(g::GrassmannianRing) = g.k * (g.n - g.k)

function deg(g::GrassmannianRing)::Integer
    if dim(g) > 20
        return _big_deg(g)
    end

    prod = factorial(g.k*(g.n-g.k))
    for i in 1:g.k
        prod *= factorial(i-1) / factorial(g.n-g.k+i-1)
    end
    return prod
end

@inline function _collect_factorial_factors!(terms::Dict{Integer, Integer}, n::Integer; func=+)
    for i in 2:n
        for (p, a) in pfactor(i)
            terms[p] = haskey(terms, p) ? func(terms[p], a) : a
        end
    end
end

function _big_deg(g::GrassmannianRing)::Integer

    terms=Dict{Integer, Integer}()
    _collect_factorial_factors!(terms, g.k * (g.n - g.k))
    for i in 2:g.k
        _collect_factorial_factors!(terms, i-1)
        _collect_factorial_factors!(terms, g.n-g.k+i-1; func=-)
    end
    _collect_factorial_factors!(terms, g.n-g.k; func=-)
    
    return prod(p^a for (p, a) in terms)
end


#=
    Pretty Printing
=#
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
