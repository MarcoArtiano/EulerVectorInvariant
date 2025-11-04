using Trixi
using Trixi: AbstractCompressibleEulerEquations, @muladd, norm
import Trixi:
    varnames,
    cons2cons,
    cons2prim,
    cons2entropy,
    entropy,
    FluxLMARS,
    boundary_condition_slip_wall,
    flux,
    max_abs_speeds,
    max_abs_speed,
    max_abs_speed_naive,
    have_nonconservative_terms,
    True, prim2cons
# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

# TODO: needs to be changed
@doc raw"""
	CompressibleEulerVectorInvariantEquations3D(gamma)
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

struct CompressibleEulerVectorInvariantEquations3D{RealT<:Real} <:
       AbstractCompressibleEulerEquations{3,6}
    p_0::RealT # reference pressure in Pa
    c_p::RealT # specific heat at constant pressure in J/(kg K)
    c_v::RealT # specific heat at constant volume in J/(kg K)
    g::RealT # gravitational acceleration in m/s²
    R::RealT # gas constant in J/(kg K)
    gamma::RealT # ratio of specific heats 
    inv_gamma_minus_one::RealT # = inv(gamma - 1); can be used to write slow divisions as fast multiplications
    K::RealT # = p_0 * (R / p_0)^gamma; scaling factor between pressure and weighted potential temperature
    stolarsky_factor::RealT # = (gamma - 1) / gamma; used in the stolarsky mean
    function CompressibleEulerVectorInvariantEquations3D(; c_p, c_v, gravity)
        c_p, c_v, g = promote(c_p, c_v, gravity)
        p_0 = 100_000
        R = c_p - c_v
        gamma = c_p / c_v
        inv_gamma_minus_one = inv(gamma - 1)
        K = p_0 * (R / p_0)^gamma
        stolarsky_factor = (gamma - 1) / gamma
        return new{typeof(c_p)}(
            p_0,
            c_p,
            c_v,
            g,
            R,
            gamma,
            inv_gamma_minus_one,
            K,
            stolarsky_factor,
        )
    end
end

function varnames(::typeof(cons2cons), ::CompressibleEulerVectorInvariantEquations3D)
    ("rho", "v1", "v2", "v3", "rho_theta", "phi")
end

have_nonconservative_terms(::CompressibleEulerVectorInvariantEquations3D) = True()

varnames(::typeof(cons2prim), ::CompressibleEulerVectorInvariantEquations3D) =
    ("rho", "v1", "v2", "v3", "p", "phi")

@inline function Trixi.boundary_condition_slip_wall(
    u_inner,
    normal_direction::AbstractVector,
    x,
    t,
    surface_flux_functions,
    equations::CompressibleEulerVectorInvariantEquations3D,
)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)

    # compute the normal velocity
		u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3] + normal[3] * u_inner[4]

    # create the "external" boundary solution state
    u_boundary = SVector(
        u_inner[1],
        u_inner[2] - 2 * u_normal * normal[1],
        u_inner[3] - 2 * u_normal * normal[2],
        u_inner[4] - 2 * u_normal * normal[3],
        u_inner[5],
	u_inner[6],
    )

    #   calculate the boundary flux
    flux, _ = surface_flux_function(u_inner, u_boundary, normal_direction, equations)
    return flux
    #= flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)
        noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
                                                     equations)
        return flux, noncons_flux =#
end

@inline function flux_energy_stable(
    u_ll,
    u_rr,
    normal_direction::AbstractVector,
    equations::CompressibleEulerVectorInvariantEquations3D)

    rho_ll, v1_ll, v2_ll, v3_ll, rho_theta_ll = u_ll
    rho_rr, v1_rr, v2_rr, v3_rr, rho_theta_rr = u_rr
    _, _, _, _, exner_ll = cons2primexner(u_ll, equations)
    _, _, _, _, exner_rr = cons2primexner(u_rr, equations)
    theta_ll = rho_theta_ll / rho_ll
    theta_rr = rho_theta_rr / rho_rr

    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] + v3_ll * normal_direction[3]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] + v3_rr * normal_direction[3]

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    kin_avg = 0.5f0 * (v1_rr * v1_rr + v2_rr * v2_rr + v3_rr * v3_rr + v1_ll * v1_ll + v2_ll * v2_ll + v3_ll * v3_ll)
    theta_avg = 0.5f0 * (theta_ll + theta_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)

    jump_v1 = v1_rr - v1_ll
    jump_v2 = v2_rr - v2_ll
    jump_v3 = v3_rr - v3_ll

    rho_v_ll = v1_ll * rho_ll * normal_direction[1] + v2_ll * rho_ll * normal_direction[2] + v3_ll * rho_ll * normal_direction[3]
    rho_v_rr = v1_rr * rho_rr * normal_direction[1] + v2_rr * rho_rr * normal_direction[2] + v3_rr * rho_rr * normal_direction[3]

    T_ll = theta_ll * exner_ll
    T_rr = theta_rr * exner_rr

    c_ll = sqrt(equations.gamma * equations.R * T_ll)
    c_rr = sqrt(equations.gamma * equations.R * T_rr)

    p_ll = exner_ll * rho_theta_ll * equations.R
    p_rr = exner_rr * rho_theta_rr * equations.R
    norm_ = norm(normal_direction)
    c = 0.5f0 * (c_ll + c_rr)
    c = 340.0
    disst = - 1 / (2 * c) * (p_rr - p_ll) * norm_
    c_adv = 0.5f0 * abs((v_dot_n_ll + v_dot_n_rr)) / norm(normal_direction) * 0 
    diss1 = 0.5f0 * c / rho_avg * (rho_v_rr - rho_v_ll) * normal_direction[1] / norm_
    diss2 = 0.5f0 * c / rho_avg * (rho_v_rr - rho_v_ll) * normal_direction[2] / norm_
    diss3 = 0.5f0 * c / rho_avg * (rho_v_rr - rho_v_ll) * normal_direction[3] / norm_
    f1 = rho_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr) + disst
    f2 = kin_avg * 0.5f0 * normal_direction[1] - diss1 -0.5f0 * c_adv / rho_avg * (rho_rr * v1_rr - rho_ll * v1_ll) *norm(normal_direction) + v1_avg/rho_avg * disst
    f3 = kin_avg * 0.5f0 * normal_direction[2] - diss2 -0.5f0 * c_adv / rho_avg * (rho_rr * v2_rr - rho_ll * v2_ll) *norm(normal_direction) + v2_avg/rho_avg * disst
    f4 = kin_avg * 0.5f0 * normal_direction[3] - diss3 -0.5f0 * c_adv / rho_avg * (rho_rr * v3_rr - rho_ll * v3_ll) *norm(normal_direction) + v3_avg/rho_avg * disst
    if f1 >= 0
        f5 = f1 * theta_ll
    else
        f5 = f1 * theta_rr
    end
    theta_grad_exner =  equations.c_p * theta_avg * (exner_rr - exner_ll)  # inv_ln_mean(1/theta_ll, 1/theta_rr)    
    vorticity_x = v2_avg * (jump_v1 * normal_direction[2] - jump_v2 * normal_direction[1]) + v3_avg * (jump_v1 * normal_direction[3] - jump_v3 * normal_direction[1])
    vorticity_y = v1_avg * (jump_v2 * normal_direction[1] - jump_v1 * normal_direction[2]) + v3_avg * (jump_v2 * normal_direction[3] - jump_v3 * normal_direction[2])
    vorticity_z = v1_avg * (jump_v3 * normal_direction[1] - jump_v1 * normal_direction[3]) + v2_avg * (jump_v3 * normal_direction[2] - jump_v2 * normal_direction[3]) 
 
    g2 = vorticity_x + theta_grad_exner * normal_direction[1]
    g3 = vorticity_y + theta_grad_exner * normal_direction[2]
    g4 = vorticity_z + theta_grad_exner * normal_direction[3]
    return SVector(f1, f2 + 0.5f0 * g2, f3 + 0.5f0 * g3, f4 + 0.5f0 * g4, f5, zero(eltype(u_ll))), SVector(f1, f2 - 0.5f0 * g2, f3 - 0.5f0 * g3, f4 - 0.5f0 * g4, f5, zero(eltype(u_ll)))
end

@inline function flux_energy_stable_mod(
    u_ll,
    u_rr,
    normal_direction::AbstractVector,
    equations::CompressibleEulerVectorInvariantEquations3D)

    rho_ll, v1_ll, v2_ll, v3_ll, rho_theta_ll = u_ll
    rho_rr, v1_rr, v2_rr, v3_rr, rho_theta_rr = u_rr
    _, _, _, _, exner_ll = cons2primexner(u_ll, equations)
    _, _, _, _, exner_rr = cons2primexner(u_rr, equations)
    theta_ll = rho_theta_ll / rho_ll
    theta_rr = rho_theta_rr / rho_rr

    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] + v3_ll * normal_direction[3]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] + v3_rr * normal_direction[3]

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    theta_avg = 0.5f0 * (theta_ll + theta_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)

    jump_v1 = v1_rr - v1_ll
    jump_v2 = v2_rr - v2_ll
    jump_v3 = v3_rr - v3_ll

    rho_v_ll = v1_ll * rho_ll * normal_direction[1] + v2_ll * rho_ll * normal_direction[2] + v3_ll * rho_ll * normal_direction[3]
    rho_v_rr = v1_rr * rho_rr * normal_direction[1] + v2_rr * rho_rr * normal_direction[2] + v3_rr * rho_rr * normal_direction[3]

    T_ll = theta_ll * exner_ll
    T_rr = theta_rr * exner_rr

    c_ll = sqrt(equations.gamma * equations.R * T_ll)
    c_rr = sqrt(equations.gamma * equations.R * T_rr)

    p_ll = exner_ll * rho_theta_ll * equations.R
    p_rr = exner_rr * rho_theta_rr * equations.R
    norm_ = norm(normal_direction)
    c = 0.5f0 * (c_ll + c_rr)
    c = 340.0
    disst = - 1 / (2 * c) * (p_rr - p_ll) * norm_
    c_adv = 0.5f0 * abs((v_dot_n_ll + v_dot_n_rr)) / norm(normal_direction)  
    diss1 = 0.5f0 * c / rho_avg * (rho_v_rr - rho_v_ll) * normal_direction[1] / norm_
    diss2 = 0.5f0 * c / rho_avg * (rho_v_rr - rho_v_ll) * normal_direction[2] / norm_
    diss3 = 0.5f0 * c / rho_avg * (rho_v_rr - rho_v_ll) * normal_direction[3] / norm_
    v_interface = 0.5f0 * (v_dot_n_ll + v_dot_n_rr) - 1 / (2 * c * rho) * (p_rr - p_ll) * norm_

    if (v_interface > 0)
	f1, f2, f3, f4, f5 = u_ll * v_interface
    else
	f1, f2, f3, f4, f5 = u_rr * v_interface
    end
    #f1 = rho_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr) + disst 
    f2 = - diss1 -0.5f0 * c_adv / rho_avg * (rho_rr * v1_rr - rho_ll * v1_ll) *norm(normal_direction) + v1_avg * disst
    f3 = - diss2 -0.5f0 * c_adv / rho_avg * (rho_rr * v2_rr - rho_ll * v2_ll) *norm(normal_direction) + v2_avg * disst
    f4 = - diss3 -0.5f0 * c_adv / rho_avg * (rho_rr * v3_rr - rho_ll * v3_ll) *norm(normal_direction) + v3_avg * disst
   # if f1 >= 0
   #     f5 = f1 * theta_ll
   # else
   #     f5 = f1 * theta_rr
   # end
    theta_grad_exner =  equations.c_p * theta_avg * (exner_rr - exner_ll)  # inv_ln_mean(1/theta_ll, 1/theta_rr)    
	
    advection = v1_avg * normal_direction[1] + v2_avg * normal_direction[2] + v3_avg * normal_direction[3] 
    advection_x = advection * jump_v1
    advection_y = advection * jump_v2
    advection_z = advection * jump_v3
    g2 = advection_x + theta_grad_exner * normal_direction[1]
    g3 = advection_y + theta_grad_exner * normal_direction[2]
    g4 = advection_z + theta_grad_exner * normal_direction[3]
    return SVector(f1, f2 + 0.5f0 * g2, f3 + 0.5f0 * g3, f4 + 0.5f0 * g4, f5, zero(eltype(u_ll))), SVector(f1, f2 - 0.5f0 * g2, f3 - 0.5f0 * g3, f4 - 0.5f0 * g4, f5, zero(eltype(u_ll)))
end

@inline function flux_lmars_mod(
    u_ll,
    u_rr,
    normal_direction::AbstractVector,
    equations::CompressibleEulerVectorInvariantEquations3D,
)

    rho_ll, v1_ll, v2_ll, v3_ll, rho_theta_ll = u_ll
    rho_rr, v1_rr, v2_rr, v3_rr, rho_theta_rr = u_rr
    _, _, _, _, exner_ll = cons2primexner(u_ll, equations)
    _, _, _, _, exner_rr = cons2primexner(u_rr, equations)
    theta_ll = rho_theta_ll / rho_ll
    theta_rr = rho_theta_rr / rho_rr
    _, _, _, _, p_ll = cons2prim(u_ll, equations)
    _, _, _, _, p_rr = cons2prim(u_rr, equations)

    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] + v3_ll * normal_direction[3]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] + v3_rr * normal_direction[3]

    v_interface = 0.5f0 * (v_dot_n_ll + v_dot_n_rr) - 1 / (2 * c * rho) * (p_rr - p_ll) * norm_

    if (v_interface > 0)
	f1, f2, f3, f4, f5 = u_ll * v_interface
    else
	f1, f2, f3, f4, f5 = u_rr * v_interface
    end

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    kin_avg = 0.5f0 * (v1_rr * v1_rr + v2_rr * v2_rr + v3_rr * v3_rr + v1_ll * v1_ll + v2_ll * v2_ll + v3_ll * v3_ll)
    theta_avg = 0.5f0 * (theta_ll + theta_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)

    jump_v1 = v1_rr - v1_ll
    jump_v2 = v2_rr - v2_ll
    jump_v3 = v3_rr - v3_ll

    rho_v_ll = v1_ll * rho_ll * normal_direction[1] + v2_ll * rho_ll * normal_direction[2] + v3_ll * rho_ll * normal_direction[3]
    rho_v_rr = v1_rr * rho_rr * normal_direction[1] + v2_rr * rho_rr * normal_direction[2] + v3_rr * rho_rr * normal_direction[3]

    T_ll = theta_ll * exner_ll
    T_rr = theta_rr * exner_rr

    c_ll = sqrt(equations.gamma * equations.R * T_ll)
    c_rr = sqrt(equations.gamma * equations.R * T_rr)

    c = 0.5f0 * (c_ll + c_rr)
    c_adv = 0.5f0 * abs((v_dot_n_ll + v_dot_n_rr)) / norm(normal_direction)
    diss1 = 0.5f0 * c / rho_avg * (rho_v_rr - rho_v_ll) * normal_direction[1] / norm(normal_direction)
    diss2 = 0.5f0 * c / rho_avg * (rho_v_rr - rho_v_ll) * normal_direction[2] / norm(normal_direction)
    diss3 = 0.5f0 * c / rho_avg * (rho_v_rr - rho_v_ll) * normal_direction[3] / norm(normal_direction)
    f1 = rho_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = kin_avg * 0.5f0 * normal_direction[1] - diss1 -0.5f0 * c_adv / rho_avg * (rho_rr * v1_rr - rho_ll * v1_ll) *norm(normal_direction)
    f3 = kin_avg * 0.5f0 * normal_direction[2] - diss2 -0.5f0 * c_adv / rho_avg * (rho_rr * v2_rr - rho_ll * v2_ll) *norm(normal_direction)
    f4 = kin_avg * 0.5f0 * normal_direction[3] - diss3 -0.5f0 * c_adv / rho_avg * (rho_rr * v3_rr - rho_ll * v3_ll) *norm(normal_direction)

    if f1 >= 0
        f5 = f1 * theta_ll
			theta_avg = theta_ll
    else
        f5 = f1 * theta_rr
			theta_avg = theta_rr
    end
    theta_grad_exner =  equations.c_p * theta_avg * (exner_rr - exner_ll)    
    vorticity_x = v2_avg * (jump_v1 * normal_direction[2] - jump_v2 * normal_direction[1]) + v3_avg * (jump_v1 * normal_direction[3] - jump_v3 * normal_direction[1])
    vorticity_y = v1_avg * (jump_v2 * normal_direction[1] - jump_v1 * normal_direction[2]) + v3_avg * (jump_v2 * normal_direction[3] - jump_v3 * normal_direction[2])
    vorticity_z = v1_avg * (jump_v3 * normal_direction[1] - jump_v1 * normal_direction[3]) + v2_avg * (jump_v3 * normal_direction[2] - jump_v2 * normal_direction[3]) 
 
    g2 = vorticity_x + theta_grad_exner * normal_direction[1]
    g3 = vorticity_y + theta_grad_exner * normal_direction[2]
    g4 = vorticity_z + theta_grad_exner * normal_direction[3]
    return SVector(f1, f2 + 0.5f0 * g2, f3 + 0.5f0 * g3, f4 + 0.5f0 * g4, f5, zero(eltype(u_ll))), SVector(f1, f2 - 0.5f0 * g2, f3 - 0.5f0 * g3, f4 - 0.5f0 * g4, f5, zero(eltype(u_ll)))
end

@inline function max_abs_speed(
    u_ll,
    u_rr,
    orientation::Integer,
    equations::CompressibleEulerVectorInvariantEquations3D,
)
    rho_ll, v1_ll, v2_ll, v3_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, v3_rr, p_rr = cons2prim(u_rr, equations)

    # Get the velocity value in the appropriate direction
    if orientation == 1
        v_ll = v1_ll
        v_rr = v1_rr
    elseif orientation == 2
        v_ll = v2_ll
        v_rr = v2_rr
    else
	v_ll = v3_ll
	v_rr = v3_rr
    end
    # Calculate sound speeds
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    return max(abs(v_ll) + c_ll, abs(v_rr) + c_rr)
end

# Convert conservative variables to primitive
@inline function cons2prim(u, equations::CompressibleEulerVectorInvariantEquations3D)
    rho, v1, v2, v3, rho_theta, phi = u

    p = equations.K * rho_theta^equations.gamma

    return SVector(rho, v1, v2, v3, p, phi)
end

# Convert primitive to conservative variables
@inline function prim2cons(prim, equations::CompressibleEulerVectorInvariantEquations3D)
    rho, v1, v2, v3, p, phi = prim
rho_theta = (p / equations.p_0)^(1 / equations.gamma) * equations.p_0 / equations.R
    return SVector(rho, v1, v2, v3, rho_theta, phi)
end

@inline function density(u, equations::CompressibleEulerVectorInvariantEquations3D)
    rho = u[1]
    return rho
end

@inline function velocity(u, equations::CompressibleEulerVectorInvariantEquations3D)
    v1 = u[2]
    v2 = u[3]
    v3 = u[4]
    return SVector(v1, v2, v3)
end

@inline function pressure(u, equations::CompressibleEulerVectorInvariantEquations3D)
    rho, v1, v2, v3, rho_theta = u
    p = equations.K * rho_theta^equations.gamma
    return p
end

@inline function cons2primexner(
    u,
    equations::CompressibleEulerVectorInvariantEquations3D,
)

    rho, v1, v2, v3, rho_theta, phi = u

    exner = (rho_theta * equations.R / equations.p_0)^(equations.R / equations.c_v)
    return SVector(rho, v1, v2, v3, exner, phi)
end

@inline function exner_pressure(
    u,
    equations::CompressibleEulerVectorInvariantEquations3D,
)

    _, _, _, _, rho_theta = u

    exner = (rho_theta * equations.R / equations.p_0)^(equations.R / equations.c_v)
    return exner
end

@inline function cons2entropy(u, equations::CompressibleEulerVectorInvariantEquations3D)
    rho, v1, v2, v3, rho_theta = u

    return SVector(0, 0, 0, 0, 0, 0)
end

end # @muladd

