using Trixi
using Trixi: AbstractCompressibleEulerEquations, @muladd, norm
import Trixi: varnames, cons2cons, cons2prim, cons2entropy, entropy, FluxLMARS, boundary_condition_slip_wall, flux, max_abs_speeds, max_abs_speed, max_abs_speed_naive, have_nonconservative_terms, True
# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# TODO: needs to be changed
@doc raw"""
	CompressibleEulerVectorInvariantEquations2D(gamma)
The compressible Euler equations
```math
\frac{\partial}{\partial t}
\begin{pmatrix}
\rho \\ \rho v_1 \\ \rho v_2 \\ \rho e
\end{pmatrix}
+
\frac{\partial}{\partial x}
\begin{pmatrix}
 \rho v_1 \\ \rho v_1^2 + p \\ \rho v_1 v_2 \\ (\rho e +p) v_1
\end{pmatrix}
+
\frac{\partial}{\partial y}
\begin{pmatrix}
\rho v_2 \\ \rho v_1 v_2 \\ \rho v_2^2 + p \\ (\rho e +p) v_2
\end{pmatrix}
=
\begin{pmatrix}
0 \\ 0 \\ 0 \\ 0
\end{pmatrix}
```
for an ideal gas with ratio of specific heats `gamma`
in two space dimensions.
Here, ``\rho`` is the density, ``v_1``, ``v_2`` the velocities, ``e`` the specific total energy **rather than** specific internal energy, and
```math
p = (\gamma - 1) \left( \rho e - \frac{1}{2} \rho (v_1^2+v_2^2) \right)
```
the pressure.
"""

struct CompressibleEulerVectorInvariantEquations2D{RealT <: Real} <:
	   AbstractCompressibleEulerEquations{2, 5}
	p_0::RealT
	c_p::RealT
	c_v::RealT
	g::RealT
	R::RealT
	gamma::RealT
	inv_gamma_minus_one::RealT
	K::RealT
	stolarsky_factor::RealT
	kappa::RealT
end

function CompressibleEulerVectorInvariantEquations2D(; g = 9.81, RealT = Float64)
	p_0 = RealT(100_000.0)
	c_p = RealT(1004.0)
	c_v = RealT(717.0)
	R = c_p - c_v
	gamma = c_p / c_v
	inv_gamma_minus_one = inv(gamma - 1.0)
	K = p_0 * (R / p_0)^gamma
	stolarsky_factor = (gamma - 1.0) / gamma
	kappa = R / c_p
	return CompressibleEulerVectorInvariantEquations2D{RealT}(p_0, c_p, c_v, g, R,
		gamma, inv_gamma_minus_one, K, stolarsky_factor, kappa)
end

function varnames(::typeof(cons2cons), ::CompressibleEulerVectorInvariantEquations2D)
	("rho", "v1", "v2", "rho_theta", "phi")
end

have_nonconservative_terms(::CompressibleEulerVectorInvariantEquations2D) = True()

varnames(::typeof(cons2prim), ::CompressibleEulerVectorInvariantEquations2D) = ("rho", "v1", "v2", "p", "phi")

@inline function source_terms_gravity(u, x, t, equations::CompressibleEulerVectorInvariantEquations2D)
	return SVector(zero(eltype(u)), zero(eltype(u)), -equations.g, zero(eltype(u)), zero(eltype(u)))
end

@inline function Trixi.boundary_condition_slip_wall(u_inner,
	normal_direction::AbstractVector,
	x, t,
	surface_flux_functions,
	equations::CompressibleEulerVectorInvariantEquations2D)
	surface_flux_function, nonconservative_flux_function = surface_flux_functions

	# normalize the outward pointing direction
	normal = normal_direction / norm(normal_direction)

	# compute the normal velocity
	u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3]

	# create the "external" boundary solution state
	u_boundary = SVector(u_inner[1],
		u_inner[2] - 2 * u_normal * normal[1],
		u_inner[3] - 2 * u_normal * normal[2],
		u_inner[4], u_inner[5])

	# calculate the boundary flux
	flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)
	noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
		equations)
	return flux, noncons_flux
end

@inline function Trixi.boundary_condition_slip_wall(u_inner, orientation,
	direction, x, t,
	surface_flux_functions,
	equations::CompressibleEulerVectorInvariantEquations2D)
	surface_flux_function, nonconservative_flux_function = surface_flux_functions

	## get the appropriate normal vector from the orientation
	if orientation == 1
		u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4], u_inner[5])
	else # orientation == 2
		u_boundary = SVector(u_inner[1], u_inner[2], -u_inner[3], u_inner[4], u_inner[5])
	end

	# Calculate boundary flux
	if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
		flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
		noncons_flux = nonconservative_flux_function(u_inner, u_boundary, orientation,
			equations)
	else # u_boundary is "left" of boundary, u_inner is "right" of boundary
		flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
		noncons_flux = nonconservative_flux_function(u_boundary, u_inner, orientation,
			equations)
	end

	return flux, noncons_flux
end

@inline function flux_surface_cons(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerVectorInvariantEquations2D)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, rho_theta_ll = u_ll
	rho_rr, v1_rr, v2_rr, rho_theta_rr = u_rr
	# rho_ll, v1_ll, v2_ll, exner_ll = cons2primexner(u_ll, equations)
	# rho_rr, v1_rr, v2_rr, exner_rr = cons2primexner(u_rr, equations)
	theta_ll = rho_theta_ll / rho_ll
	theta_rr = rho_theta_rr / rho_rr

	# Average each factor of products in flux
	rho_avg = 0.5f0 * (rho_ll + rho_rr)
	# v1_avg = 0.5f0 * (v1_ll + v1_rr)
	# v2_avg = 0.5f0 * (v2_ll + v2_rr)
	# exner_avg = 0.5f0 * (exner_ll + exner_rr)
	theta_avg = 0.5f0 * (theta_ll + theta_rr)
	kin_avg = 0.5f0 * (v1_rr * v1_rr + v2_rr * v2_rr + v1_ll * v1_ll + v2_ll * v2_ll)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

	## According to Kieran notes I should use the average of the momentum in the density and potential temperature fluxes?
	f1 = rho_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = kin_avg * 0.5f0 * normal_direction[1]
	f3 = kin_avg * 0.5f0 * normal_direction[2]
	f4 = f1 * theta_avg

	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux_surface_noncons(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerVectorInvariantEquations2D)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, rho_theta_ll = u_ll
	rho_rr, v1_rr, v2_rr, rho_theta_rr = u_rr
	rho_ll, v1_ll, v2_ll, exner_ll = cons2primexner(u_ll, equations)
	rho_rr, v1_rr, v2_rr, exner_rr = cons2primexner(u_rr, equations)
	theta_ll = rho_theta_ll / rho_ll
	# theta_rr = rho_theta_rr / rho_rr

	# Average each factor of products in flux
	# rho_avg = 0.5f0 * (rho_ll + rho_rr)
	# v1_avg = 0.5f0 * (v1_ll + v1_rr)
	# v2_avg = 0.5f0 * (v2_ll + v2_rr)
	# exner_avg = 0.5f0 * (exner_ll + exner_rr)
	# theta_avg = 0.5f0 * (theta_ll + theta_rr)
	# kin_avg = 0.5f0 * (v1_rr * v1_rr + v2_rr * v2_rr + v1_ll * v1_ll + v2_ll * v2_ll)
	# v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	# v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

	jump_v1 = v1_rr - v1_ll
	jump_v2 = v2_rr - v1_ll

	f1 = 0.0
	f2 = v2_ll * jump_v1 * normal_direction[2] - v2_ll * jump_v2 * normal_direction[1] + equations.c_p * theta_ll * (exner_rr - exner_ll) * normal_direction[1]
	f3 = v1_ll * jump_v2 * normal_direction[1] - v1_ll * jump_v1 * normal_direction[2] + equations.c_p * theta_ll * (exner_rr - exner_ll) * normal_direction[2]
	f4 = 0.0
	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux_surface_cons_upwind(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerVectorInvariantEquations2D)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, rho_theta_ll = u_ll
	rho_rr, v1_rr, v2_rr, rho_theta_rr = u_rr
	theta_ll = rho_theta_ll / rho_ll
	theta_rr = rho_theta_rr / rho_rr

	# Average each factor of products in flux
	rho_avg = 0.5f0 * (rho_ll + rho_rr)

	kin_avg = 0.5f0 * (v1_rr * v1_rr + v2_rr * v2_rr + v1_ll * v1_ll + v2_ll * v2_ll)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

	rho_v_rr = v1_ll * rho_ll * normal_direction[1] + v2_ll * rho_ll * normal_direction[2]
	rho_v_ll = v1_rr * rho_rr * normal_direction[1] + v2_rr * rho_rr * normal_direction[2]
	RealT = eltype(u_ll)
	c = RealT(340.0)

	c_adv = 0.5f0 * abs((v_dot_n_ll + v_dot_n_rr))/norm(normal_direction)
	diss1 = c / (2 * rho_avg) * (rho_v_rr - rho_v_ll) * normal_direction[1] / norm(normal_direction)^2
	diss2 = c / (2 * rho_avg) * (rho_v_rr - rho_v_ll) * normal_direction[2] / norm(normal_direction)^2
	## According to Kieran notes I should use the average of the momentum in the density and potential temperature fluxes?
	f1 = rho_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = kin_avg * 0.5f0 * normal_direction[1] - diss1 - 0.5f0 * c_adv/rho_avg * (rho_rr * v1_rr - rho_ll * v1_ll) * norm(normal_direction)
	f3 = kin_avg * 0.5f0 * normal_direction[2] - diss2 - 0.5f0 * c_adv/rho_avg * (rho_rr * v2_rr - rho_ll * v2_ll) * norm(normal_direction)

	if f1 >= 0
		f4 = f1 * theta_ll
	else
		f4 = f1 * theta_rr
	end

	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux_surface_noncons_upwind(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerVectorInvariantEquations2D)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, rho_theta_ll = u_ll
	rho_rr, v1_rr, v2_rr, rho_theta_rr = u_rr
	rho_ll, v1_ll, v2_ll, exner_ll = cons2primexner(u_ll, equations)
	rho_rr, v1_rr, v2_rr, exner_rr = cons2primexner(u_rr, equations)
	theta_ll = rho_theta_ll / rho_ll
	theta_rr = rho_theta_rr / rho_rr

	# Average each factor of products in flux
	rho_avg = 0.5f0 * (rho_ll + rho_rr)
	# v1_avg = 0.5f0 * (v1_ll + v1_rr)
	# v2_avg = 0.5f0 * (v2_ll + v2_rr)
	# exner_avg = 0.5f0 * (exner_ll + exner_rr)
	# theta_avg = 0.5f0 * (theta_ll + theta_rr)
	# kin_avg = 0.5f0 * (v1_rr * v1_rr + v2_rr * v2_rr + v1_ll * v1_ll + v2_ll * v2_ll)
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

	f1 = 0.0
	f2 = v2_ll * jump_v1 * normal_direction[2] - v2_ll * jump_v2 * normal_direction[1] + equations.c_p * theta * (exner_rr - exner_ll) * normal_direction[1]
	f3 = v1_ll * jump_v2 * normal_direction[1] - v1_ll * jump_v1 * normal_direction[2] + equations.c_p * theta * (exner_rr - exner_ll) * normal_direction[2]
	f4 = 0.0
	return SVector(f1, f2, f3, f4, 0)

end

@inline function flux_volume_cons(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerVectorInvariantEquations2D)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, rho_theta_ll = u_ll
	rho_rr, v1_rr, v2_rr, rho_theta_rr = u_rr
	theta_ll = rho_theta_ll / rho_ll
	theta_rr = rho_theta_rr / rho_rr

	# Average each factor of products in flux
	rho_avg = 0.5f0 * (rho_ll + rho_rr)
	theta_avg = 0.5f0 * (theta_ll + theta_rr)
	kin_avg = 0.5f0 * (v1_rr * v1_rr + v2_rr * v2_rr + v1_ll * v1_ll + v2_ll * v2_ll)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

	## According to Kieran notes I should use the average of the momentum in the density and potential temperature fluxes?
	f1 = rho_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = kin_avg * 0.5f0 * normal_direction[1]
	f3 = kin_avg * 0.5f0 * normal_direction[2]
	f4 = f1 * theta_avg

	return SVector(f1, f2, f3, f4, 0)
end


@inline function flux_volume_noncons(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerVectorInvariantEquations2D)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, rho_theta_ll, phi_ll = u_ll
	rho_rr, v1_rr, v2_rr, rho_theta_rr, phi_rr = u_rr
	rho_ll, v1_ll, v2_ll, exner_ll = cons2primexner(u_ll, equations)
	rho_rr, v1_rr, v2_rr, exner_rr = cons2primexner(u_rr, equations)
	theta_ll = rho_theta_ll/rho_ll
	# Average each factor of products in flux

	jump_v1 = v1_rr - v1_ll
	jump_v2 = v2_rr - v1_ll
	phi_jump = phi_rr - phi_ll
	f1 = 0.0
	f2 = v2_ll * jump_v1 * normal_direction[2] - v2_ll * jump_v2 * normal_direction[1] + equations.c_p * theta_ll * (exner_rr - exner_ll) * normal_direction[1] + equations.g * phi_jump * normal_direction[1]
	f3 = v1_ll * jump_v2 * normal_direction[1] - v1_ll * jump_v1 * normal_direction[2] + equations.c_p * theta_ll * (exner_rr - exner_ll) * normal_direction[2] + equations.g * phi_jump * normal_direction[2]
	#	f4 = theta_ll * (rho_rr * v1_rr - rho_ll * v1_ll) * normal_direction[1] * 0.5 + rho_ll * v1_ll * (theta_rr - theta_ll) * normal_direction[1] *0.5 +  theta_ll * (rho_rr * v2_rr - rho_ll * v2_ll) * normal_direction[2] * 0.5 + rho_ll * v2_ll * (theta_rr - theta_ll) * normal_direction[2] * 0.5  
	f4 = 0.0
	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux_zero(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerVectorInvariantEquations2D)

	return SVector(0, 0, 0, 0, 0)
end



@inline function max_abs_speed(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerVectorInvariantEquations2D)
	rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

	# Get the velocity value in the appropriate direction
	if orientation == 1
		v_ll = v1_ll
		v_rr = v1_rr
	else # orientation == 2
		v_ll = v2_ll
		v_rr = v2_rr
	end
	# Calculate sound speeds
	c_ll = sqrt(equations.gamma * p_ll / rho_ll)
	c_rr = sqrt(equations.gamma * p_rr / rho_rr)

	return max(abs(v_ll) + c_ll, abs(v_rr) + c_rr)
end

@inline function min_max_speed_naive(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerVectorInvariantEquations2D)
	rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

	if orientation == 1 # x-direction
		λ_min = v1_ll - sqrt(equations.gamma * p_ll / rho_ll)
		λ_max = v1_rr + sqrt(equations.gamma * p_rr / rho_rr)
	else # y-direction
		λ_min = v2_ll - sqrt(equations.gamma * p_ll / rho_ll)
		λ_max = v2_rr + sqrt(equations.gamma * p_rr / rho_rr)
	end

	return λ_min, λ_max
end

@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerVectorInvariantEquations2D)
	rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

	# Calculate normal velocities and sound speed
	# left
	v_ll = (v1_ll * normal_direction[1]
			+
			v2_ll * normal_direction[2])
	c_ll = sqrt(equations.gamma * p_ll / rho_ll)
	# right
	v_rr = (v1_rr * normal_direction[1]
			+
			v2_rr * normal_direction[2])
	c_rr = sqrt(equations.gamma * p_rr / rho_rr)

	return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction)
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function max_abs_speed(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerVectorInvariantEquations2D)
	rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

	# Calculate normal velocities and sound speeds
	# left
	v_ll = (v1_ll * normal_direction[1]
			+
			v2_ll * normal_direction[2])
	c_ll = sqrt(equations.gamma * p_ll / rho_ll)
	# right
	v_rr = (v1_rr * normal_direction[1]
			+
			v2_rr * normal_direction[2])
	c_rr = sqrt(equations.gamma * p_rr / rho_rr)

	norm_ = norm(normal_direction)
	return max(abs(v_ll) + c_ll * norm_,
		abs(v_rr) + c_rr * norm_)
end

@inline function max_abs_speeds(u, equations::CompressibleEulerVectorInvariantEquations2D)
	rho, v1, v2, p = cons2prim(u, equations)
	c = sqrt(equations.gamma * p / rho)

	return abs(v1) + c, abs(v2) + c
end

# Convert conservative variables to primitive
@inline function cons2prim(u, equations::CompressibleEulerVectorInvariantEquations2D)
	rho, v1, v2, rho_theta, phi = u

	p = equations.K * rho_theta^equations.gamma

	return SVector(rho, v1, v2, p, phi)
end

# Convert primitive to conservative variables
@inline function prim2cons(prim, equations::CompressibleEulerVectorInvariantEquations2D)
	rho, v1, v2, p, phi = prim
	rho_theta = (p / equations.p_0)^(1 / equations.gamma) * equations.p_0 / equations.R
	return SVector(rho, v1, v2, rho_theta, phi)
end

@inline function density(u, equations::CompressibleEulerVectorInvariantEquations2D)
	rho = u[1]
	return rho
end

@inline function velocity(u, equations::CompressibleEulerVectorInvariantEquations2D)
	v1 = u[2]
	v2 = u[3]
	return SVector(v1, v2)
end

@inline function pressure(u, equations::CompressibleEulerVectorInvariantEquations2D)
	rho, v1, v2, rho_theta = u
	p = equations.K * rho_theta^equations.gamma
	return p
end

@inline function cons2primexner(u, equations::CompressibleEulerVectorInvariantEquations2D)

	rho, v1, v2, rho_theta, phi = u

	exner = (rho_theta * equations.R / equations.p_0)^(equations.R / equations.c_v)
	return SVector(rho, v1, v2, exner, phi)
end

@inline function exner_pressure(u, equations::CompressibleEulerVectorInvariantEquations2D)

	_, _, _, rho_theta = u

	exner = (rho_theta * equations.R / equations.p_0)^(equations.R / equations.c_v)
	return exner
end

@inline function cons2entropy(u, equations::CompressibleEulerVectorInvariantEquations2D)
	rho, v1, v2, rho_theta = u

	return SVector(0, 0, 0, 0, 0)
end

@inline function well_balanced_v1(u, equations::CompressibleEulerVectorInvariantEquations2D)
	rho, v1, v2, rho_theta, _ = u
	return abs(v1)
end

@inline function well_balanced_v2(u, equations::CompressibleEulerVectorInvariantEquations2D)
	rho, v1, v2, rho_theta, _ = u
	return abs(v2)
end

end # @muladd
