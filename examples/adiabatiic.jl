using OrdinaryDiffEqLowStorageRK
using Invariant
using Invariant.Trixi
using Plots
using OrdinaryDiffEq
using DoubleFloats

function initial_condition_adiabatic(x, t, equations::CompressibleEulerVectorInvariantEquations2D)
    RealT = eltype(x)
    g = RealT(9.81)
    c_p = RealT(1004.0)
    c_v = RealT(717.0)
    # center of perturbation
    p0 = RealT(100_000)
    # perturbation in potential temperature
    R = c_p - c_v    # gas constant (dry air)
    potential_temperature = RealT(300.0)

    # Exner pressure, solves hydrostatic equation for x[2]
    exner = RealT(1) - g / (c_p * potential_temperature) * x[2]

    # pressure
    R = c_p - c_v    # gas constant (dry air)
    p = p0 * exner^(c_p / R)

    # temperature
    T = potential_temperature * exner
    # density
    rho = p / (R * T)
    v1 = RealT(0.0)
    v2 = RealT(0.0)

    return Invariant.prim2cons(SVector(rho, v1, v2, p, x[2]), equations)
end

RealT = Float64
T = RealT(100)
equations = CompressibleEulerVectorInvariantEquations2D()

boundary_conditions = Dict(:y_neg => boundary_condition_slip_wall,
		           		   :y_pos => boundary_condition_slip_wall)

polydeg = 2
basis = LobattoLegendreBasis(RealT, polydeg)

# surface_flux_diss = FluxPlusDissipation(flux_surface_cons,DissipationLocalLaxFriedrichs(max_abs_speed_naive))
# surface_flux = (flux_surface_cons, flux_surface_noncons)


@inline function flux_surface_noncons_wb(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerVectorInvariantEquations2D)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, rho_theta_ll = u_ll
	rho_rr, v1_rr, v2_rr, rho_theta_rr = u_rr
	rho_ll, v1_ll, v2_ll, exner_ll = Invariant.cons2primexner(u_ll, equations)
	rho_rr, v1_rr, v2_rr, exner_rr = Invariant.cons2primexner(u_rr, equations)
	theta_ll = rho_theta_ll / rho_ll
	theta_rr = rho_theta_rr / rho_rr

	# Average each factor of products in flux
	rho_avg = 0.5f0 * (rho_ll + rho_rr)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
	
	jump_v1 = v1_rr - v1_ll
	jump_v2 = v2_rr - v1_ll
	f1 = rho_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	if f1 >= 0 
		theta = theta_ll
	else
		theta = theta_rr
	end

    theta = 0.5f0 * (theta_ll + theta_rr)
	theta = Trixi.inv_ln_mean(1/theta_ll, 1/theta_rr)
	#theta = 0.0
    f1 = 0.0
	f2 = v2_ll * jump_v1 * normal_direction[2] -v2_ll * jump_v2 * normal_direction[1] + equations.c_p * theta * (exner_rr - exner_ll) * normal_direction[1]
	f3 = v1_ll * jump_v2 * normal_direction[1] -v1_ll * jump_v1 * normal_direction[2] + equations.c_p * theta * (exner_rr - exner_ll) * normal_direction[2]
	f4 = 0.0
	return SVector(f1, f2, f3, f4, 0)

end

@inline function flux_volume_wb(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerVectorInvariantEquations2D)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, rho_theta_ll, phi_ll = u_ll
	rho_rr, v1_rr, v2_rr, rho_theta_rr, phi_rr = u_rr
	rho_ll, v1_ll, v2_ll, exner_ll = Invariant.cons2primexner(u_ll, equations)
	rho_rr, v1_rr, v2_rr, exner_rr = Invariant.cons2primexner(u_rr, equations)
	theta_ll = rho_theta_ll/rho_ll
	theta_rr = rho_theta_rr/rho_rr

	# Average each factor of products in flux
	jump_v1 = v1_rr - v1_ll
	jump_v2 = v2_rr - v1_ll
	phi_jump = phi_rr - phi_ll
	theta_avg = (theta_ll + theta_rr)*0.5f0
    f1 = 0.0
	f2 = v2_ll * jump_v1 * normal_direction[2] -v2_ll * jump_v2 * normal_direction[1] + equations.c_p * theta_avg * (exner_rr - exner_ll) * normal_direction[1] + equations.g * phi_jump * normal_direction[1]
	f3 = v1_ll * jump_v2 * normal_direction[1] -v1_ll * jump_v1 * normal_direction[2] + equations.c_p * theta_avg * (exner_rr - exner_ll) * normal_direction[2] + equations.g * phi_jump * normal_direction[2]
	f4 = 0.0
	return SVector(f1, f2, f3, f4, 0)
end

surface_flux = (flux_surface_cons_upwind, flux_surface_noncons_wb)
volume_flux = (flux_volume_cons, flux_volume_wb)
surface_flux = (flux_surface_cons_upwind, flux_zero)
volume_flux = (flux_volume_cons, flux_volume_wb)
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

trees_per_dimension = (16, 16)

function mapping(xi, eta)
    x = xi + RealT(0.1) * sinpi(xi) * sinpi(eta)
    y = eta + RealT(0.1) * sinpi(xi) * sinpi(eta)
    return SVector(RealT(1000) * RealT(0.5) * (RealT(1) + x), RealT(1000) * RealT(0.5) * (RealT(1) + y))
end

mesh = P4estMesh(trees_per_dimension, polydeg = polydeg, 
mapping=mapping, periodicity = (false,false), initial_refinement_level=0, RealT = RealT)

coordinates_min = (0.0, 0.0)
coordinates_max = (1_000.0, 1_000.0)

#cells_per_dimension = (16, 16)
#mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
#	periodicity = (true, false))

trees_per_dimension = (16, 16)

mesh = P4estMesh(trees_per_dimension, polydeg = polydeg,
	coordinates_min = coordinates_min, coordinates_max = coordinates_max,
	periodicity = (false, false), initial_refinement_level = 0)

boundary_conditions =  Dict( :x_pos => boundary_condition_slip_wall, 
                             :x_neg => boundary_condition_slip_wall,
                             :y_pos => boundary_condition_slip_wall,
                             :y_neg => boundary_condition_slip_wall)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_isothermal, solver, boundary_conditions=boundary_conditions)

dt = RealT(0.01)
tspan = (zero(RealT), 2 * dt)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval=analysis_interval, extra_analysis_integrals=(Invariant.well_balanced_v1, Invariant.well_balanced_v2))#, save_analysis = true, output_directory = pwd()*"/test_cases/balance/out", analysis_filename = "well_balancing_$(RealT).dat")

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)


sol = solve(ode, Euler(); 
                 dt = dt, ode_default_options()..., callback = callbacks, adaptive = false)

summary_callback();