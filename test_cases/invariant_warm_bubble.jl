using OrdinaryDiffEqLowStorageRK
using Invariant
using Invariant.Trixi
using Plots
using OrdinaryDiffEqSSPRK

# Initial condition
function initial_condition_warm_bubble(
    x,
    t,
    equations::CompressibleEulerVectorInvariantEquations2D,
)
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
    return SVector(rho, v1, v2, rho * potential_temperature, g * x[2])
end

###############################################################################
# semidiscretization of the compressible Euler equations
equations = CompressibleEulerVectorInvariantEquations2D(c_p = 1007, c_v = 717, gravity = 9.81)

boundary_conditions =
    Dict(:y_neg => boundary_condition_slip_wall, :y_pos => boundary_condition_slip_wall)

polydeg = 3
basis = LobattoLegendreBasis(polydeg)

# surface_flux_diss = FluxPlusDissipation(flux_surface_cons,DissipationLocalLaxFriedrichs(max_abs_speed_naive))
# surface_flux = (flux_surface_cons, flux_surface_noncons)

@inline function flux_volume_wb(
    u_ll,
    u_rr,
    normal_direction::AbstractVector,
    equations::CompressibleEulerVectorInvariantEquations2D,
)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, rho_theta_ll, phi_ll = u_ll
    rho_rr, v1_rr, v2_rr, rho_theta_rr, phi_rr = u_rr
    rho_ll, v1_ll, v2_ll, exner_ll = Invariant.cons2primexner(u_ll, equations)
    rho_rr, v1_rr, v2_rr, exner_rr = Invariant.cons2primexner(u_rr, equations)
    theta_ll = rho_theta_ll / rho_ll
    theta_rr = rho_theta_rr / rho_rr

    # Average each factor of products in flux

    jump_v1 = v1_rr - v1_ll
    jump_v2 = v2_rr - v1_ll
    phi_jump = phi_rr - phi_ll
    theta_avg = (theta_ll + theta_rr) * 0.5f0
    theta_avg = Trixi.inv_ln_mean(1 / theta_ll, 1 / theta_rr)
    f1 = 0.0
    f2 =
        v2_ll * jump_v1 * normal_direction[2] - v2_ll * jump_v2 * normal_direction[1] +
        theta_avg * (exner_rr - exner_ll) * normal_direction[1] +
        equations.g * phi_jump * normal_direction[1]
    f3 =
        v1_ll * jump_v2 * normal_direction[1] - v1_ll * jump_v1 * normal_direction[2] +
        theta_avg * (exner_rr - exner_ll) * normal_direction[2] +
        equations.g * phi_jump * normal_direction[2]
    #	f4 = theta_ll * (rho_rr * v1_rr - rho_ll * v1_ll) * normal_direction[1] * 0.5 + rho_ll * v1_ll * (theta_rr - theta_ll) * normal_direction[1] *0.5 +  theta_ll * (rho_rr * v2_rr - rho_ll * v2_ll) * normal_direction[2] * 0.5 + rho_ll * v2_ll * (theta_rr - theta_ll) * normal_direction[2] * 0.5  
    f4 = 0.0
    return SVector(f1, f2, f3, f4, 0)
end

function flux_zero_n(u_ll, u_rr, normal_or_orientation, equations)
    return zero(u_ll), zero(u_rr)
end

#surface_flux = (flux_surface_cons_upwind, flux_surface_noncons_upwind)
#surface_flux = (flux_surface_total_upwind, flux_zero)
surface_flux = (flux_energy_stable, flux_zero)
volume_flux = (flux_invariant_turbo, flux_zero)
#volume_flux = (flux_volume_cons, flux_volume_noncons)
#volume_flux = (flux_zero, flux_zero)

volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (0.0, 0.0)
coordinates_max = (20_000.0, 10_000.0)

#cells_per_dimension = (16, 16)
#mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
#	periodicity = (true, false))
# 32, 16
trees_per_dimension = (32, 16)

mesh = P4estMesh(
    trees_per_dimension,
    polydeg = polydeg,
    coordinates_min = coordinates_min,
    coordinates_max = coordinates_max,
    periodicity = (true, false),
    initial_refinement_level = 0,
)

semi = SemidiscretizationHyperbolic(
    mesh,
    equations,
    initial_condition_warm_bubble,
    solver,
    boundary_conditions = boundary_conditions,
)

###############################################################################
# ODE solvers, callbacks etc.
tspan = (0.0, 1000.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000

analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)

###############################################################################
# run the simulation
sol =
    solve(ode, SSPRK43(); maxiters = 1.0e7, ode_default_options()..., callback = callbacks);
