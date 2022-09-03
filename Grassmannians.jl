using AbstractAlgebra
import Base.*
import Base.+
import Base.-
import Base.==

include("schubert.jl")

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

mymat = g[a b;c d]
det(mymat)
m = _giambelli_matrix(g, Partition([3,2]))