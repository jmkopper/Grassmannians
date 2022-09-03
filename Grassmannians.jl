using AbstractAlgebra
import Base.*
import Base.+
import Base.==

include("schubert.jl")

function _valid_pieri_summand(pieripartition::Generic.Partition, p::Generic.Partition)::Bool
    for i in 1:min(length(pieripartition), length(p))
        if pieripartition[i] > p[i] || (i > 1 && pieripartition[i-1] < p[i])
            return false
        end
    end
    return true
end

function _is_pieri(p::Generic.Partition)::Bool
  return (length(p) <= 1)
end

function pieri_prod(p::Generic.Partition, q::Generic.Partition)::Vector{Generic.Partition}
    v_n = p.n + q.n
    valid_partitions = []
    for part in Generic.partitions(v_n)
        if _valid_pieri_summand(q, part)
            push!(valid_partitions, part)
        end
    end

    return valid_partitions
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

    return pieri_prod(p, q)
end

g = GrassmannianRing(3, 6)
terms_a = Dict(Partition([3,2,1])=>1, Partition([3])=>2)
a = SchubertCycle(3, 6, terms_a, g)

terms_b = Dict(Partition([1])=>3)
b = SchubertCycle(3, 6, terms_b, g)
println(a*b)
println(zero(g))
println(g())
c = g([3,3,3])
a = g([3])
println(c)
println(c*a)

a = g([3])
b = g([2])
c = g([2])
d = g([1])

mymat = [a b;c d]
det(mymat)