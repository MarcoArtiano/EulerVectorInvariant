module Invariant

using Trixi

include("compressible_euler_vectorinvariant_2d.jl")

include("dg_invariant_2d.jl")

export CompressibleEulerVectorInvariantEquations2D

end # module Invariant
