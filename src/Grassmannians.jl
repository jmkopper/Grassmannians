module Grassmannians

using AbstractAlgebra

import Base: *, +, -, ==, ^
import AbstractAlgebra: zero, one, elem_type, parent_type, parent
import Primes: factor as pfactor

include("types.jl")
include("schubert.jl")

export GrassmannianRing, SchubertCycle
export +, -, *, ==, ^
export one, is_zero, zero
export deg, dim

end # module Grassmannians