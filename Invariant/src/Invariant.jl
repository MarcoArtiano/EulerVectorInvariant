module Invariant

using Trixi
using InteractiveUtils
include("compressible_euler_vectorinvariant_2d.jl")
include("compressible_euler_vectorinvariant_3d.jl")

# include("dg_invariant_2d.jl")
include("noncons_kernel_2d.jl")
include("noncons_kernel_3d.jl")
export CompressibleEulerVectorInvariantEquations2D
export CompressibleEulerVectorInvariantEquations3D

export flux_zero, flux_invariant_turbo, flux_energy_stable, flux_energy_stable_mod, flux_invariant_adv_turbo

export well_balanced_v1, well_balanced_v2

end # module Invariant
