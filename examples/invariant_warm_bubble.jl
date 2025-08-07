using OrdinaryDiffEqLowStorageRK
using Invariant
using Invariant.Trixi
using Plots
using OrdinaryDiffEq

# Initial condition
function initial_condition_warm_bubble(x, t, equations::CompressibleEulerVectorInvariantEquations2D)
	g = equations.g
	c_p = equations.c_p
	c_v = equations.c_v
	# center of perturbation
	center_x = 10000.0
	center_z = 2000.0
	# radius of perturbation
	radius = 2000.0
	# distance of current x to center of perturbation
	r = sqrt((x[1] - center_x)^2 + (x[2] - center_z)^2)

	# perturbation in potential temperature
	potential_temperature_ref = 300.0
	potential_temperature_perturbation = 0.0
	if r <= radius
		potential_temperature_perturbation = 2 * cospi(0.5 * r / radius)^2
	end
	potential_temperature = potential_temperature_ref + potential_temperature_perturbation

	# Exner pressure, solves hydrostatic equation for x[2]
	exner = 1 - g / (c_p * potential_temperature) * x[2]

	# pressure
	p_0 = 100_000.0  # reference pressure
	R = c_p - c_v    # gas constant (dry air)
	p = p_0 * exner^(c_p / R)

	# temperature
	T = potential_temperature * exner

	# density
	rho = p / (R * T)

	v1 = 20.0
	v2 = 0.0
	return SVector(rho, v1, v2, rho * potential_temperature)
end

###############################################################################
# semidiscretization of the compressible Euler equations
equations = CompressibleEulerVectorInvariantEquations2D()

boundary_conditions = Dict(:y_neg => boundary_condition_slip_wall,
		           		   :y_pos => boundary_condition_slip_wall)

polydeg = 3
basis = LobattoLegendreBasis(polydeg)

# surface_flux_diss = FluxPlusDissipation(flux_surface_cons,DissipationLocalLaxFriedrichs(max_abs_speed_naive))
# surface_flux = (flux_surface_cons, flux_surface_noncons)

surface_flux = (flux_surface_cons_upwind, flux_surface_noncons_upwind)
volume_flux = (flux_volume_cons, flux_volume_noncons)
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (0.0, 0.0)
coordinates_max = (20_000.0, 10_000.0)

#cells_per_dimension = (16, 16)
#mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
#	periodicity = (true, false))

trees_per_dimension = (64, 32)

mesh = P4estMesh(trees_per_dimension, polydeg = polydeg,
	coordinates_min = coordinates_min, coordinates_max = coordinates_max,
	periodicity = (true, false), initial_refinement_level = 0)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_warm_bubble, solver,
	source_terms = source_terms_gravity,
	boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.
dt = 1.0e-3
tspan = (0.0, 1000)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
	analysis_callback,
	alive_callback,
	stepsize_callback
	)

###############################################################################
# run the simulation
sol = solve(ode,
	#Euler();
	CarpenterKennedy2N54(williamson_condition = false);
	maxiters = 1.0e7,
	dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
	ode_default_options()..., callback = callbacks);

