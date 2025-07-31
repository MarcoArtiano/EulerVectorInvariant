using Trixi
using Trixi: AbstractCompressibleEulerEquations
import Trixi: varnames, cons2cons, cons2prim, cons2entropy, entropy, FluxLMARS, boundary_condition_slip_wall, flux, max_abs_speed, max_abs_speed_naive, @muladd
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

struct CompressibleEulerVectorInvariantEquations2D{RealT<:Real} <:
       AbstractCompressibleEulerEquations{2,4}
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
    p_0 = 100_000.0
    c_p = 1004.0
    c_v = 717.0
    R = c_p - c_v
    gamma = c_p / c_v
    inv_gamma_minus_one = inv(gamma - 1.0)
    K = p_0 * (R / p_0)^gamma
    stolarsky_factor = (gamma - 1.0) / gamma
    kappa = R/c_p
    return CompressibleEulerVectorInvariantEquations2D{RealT}(p_0, c_p, c_v, g, R,
        gamma, inv_gamma_minus_one, K, stolarsky_factor, kappa)
end

function varnames(::typeof(cons2cons), ::CompressibleEulerVectorInvariantEquations2D)
    ("rho", "v1", "v2", "rho_theta")
end
#TODO: maybe we could put exner...
varnames(::typeof(cons2prim), ::CompressibleEulerVectorInvariantEquations2D) = ("rho", "v1", "v2", "p")

@inline function source_terms_gravity(u, x, t, equations::CompressibleEulerVectorInvariantEquations2D)
    rho, _, _, _ = u
    return SVector(zero(eltype(u)), zero(eltype(u)), -equations.g * rho, 0.0)
end

# TODO:
# Calculate 1D flux for a single point in the normal direction
# Note, this directional vector is not normalized
@inline function flux(u, normal_direction::AbstractVector,
                      equations::CompressibleEulerVectorInvariantEquations2D)
    rho_e = last(u)
    rho, v1, v2, p = cons2prim(u, equations)

    v_normal = v1 * normal_direction[1] + v2 * normal_direction[2]
    rho_v_normal = rho * v_normal
    f1 = rho_v_normal
    f2 = rho_v_normal * v1 + p * normal_direction[1]
    f3 = rho_v_normal * v2 + p * normal_direction[2]
    f4 = (rho_e + p) * v_normal
    return SVector(f1, f2, f3, f4)
end

@inline function flux(u, orientation::Integer, equations::CompressibleEulerVectorInvariant2D)
    rho, v1, v2, rho_theta = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5f0 * (rho_v1 * v1 + rho_v2 * v2))
    if orientation == 1
        f1 = rho * v1
        f2 = 0
        f3 = 0
        f4 = rho_theta * v1
    else
        f1 = rho * v2
        f2 = 0
        f3 = 0
        f4 = rho_theta * v2
    end
    return SVector(f1, f2, f3, f4)
end


# TODO: left for reference
@inline function flux_kennedy_gruber(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::CompressibleEulerVectorInvariantEquations2D)
    # Unpack left and right state
    rho_e_ll = last(u_ll)
    rho_e_rr = last(u_rr)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v_dot_n_avg = v1_avg * normal_direction[1] + v2_avg * normal_direction[2]
    p_avg = 0.5f0 * (p_ll + p_rr)
    e_avg = 0.5f0 * (rho_e_ll / rho_ll + rho_e_rr / rho_rr)

    # Calculate fluxes depending on normal_direction
    f1 = rho_avg * v_dot_n_avg
    f2 = f1 * v1_avg + p_avg * normal_direction[1]
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * e_avg + p_avg * v_dot_n_avg

    return SVector(f1, f2, f3, f4)
end

@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector,
    x, t,
    surface_flux_function,
    equations::CompressibleEulerVectorInvariantEquations2D)
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)

    # compute the normal velocity
    u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3]

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
        u_inner[2]- 2 * u_normal * normal[1],
        u_inner[3] - 2 * u_normal * normal[2],
        u_inner[4])

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)

    return flux
end


## TODO: Left for reference
@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, normal_direction::AbstractVector,
                                         equations::CompressibleEulerVectorInvariantEquations2D)
    c = flux_lmars.speed_of_sound

    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # Note that this is the same as computing v_ll and v_rr with a normalized normal vector
    # and then multiplying v by `norm_` again, but this version is slightly faster.
    norm_ = norm(normal_direction)

    rho = 0.5f0 * (rho_ll + rho_rr)
    p = 0.5f0 * (p_ll + p_rr) - 0.5f0 * c * rho * (v_rr - v_ll) / norm_
    v = 0.5f0 * (v_ll + v_rr) - 1 / (2 * c * rho) * (p_rr - p_ll) * norm_

    # We treat the energy term analogous to the potential temperature term in the paper by
    # Chen et al., i.e. we use p_ll and p_rr, and not p
    if v >= 0
        f1, f2, f3, f4 = u_ll * v
        f4 = f4 + p_ll * v
    else
        f1, f2, f3, f4 = u_rr * v
        f4 = f4 + p_rr * v
    end

    return SVector(f1,
                   f2 + p * normal_direction[1],
                   f3 + p * normal_direction[2],
                   f4)
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
    rho, v1, v2, rho_theta = u

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    p = equations.K *rho_theta^equations.gamma

    return SVector(rho, v1, v2, p)
end

# Convert primitive to conservative variables
@inline function prim2cons(prim, equations::CompressibleEulerVectorInvariantEquations2D)
    rho, v1, v2, p = prim
    rho_v1 = rho * v1
    rho_v2 = rho * v2
    rho_theta = (p / equations.p_0)^(1 / equations.gamma) * equations.p_0 / equations.R
    return SVector(rho, rho_v1, rho_v2, rho_theta)
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
    p = equations.K *rho_theta^equations.gamma
    return p
end

@inline function cons2primexner(u, equations::CompressibleEulerVectorInvariantEquations2D)

    rho, v1, v2, rho_theta = u

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    exner = equations.c_p * (rho_theta*equations.R/equations.p_0)^(equations.R/equations.c_v)
    return SVector(rho, v1, v2, exner)
end

@inline function exner_pressure(u, equations::CompressibleEulerVectorInvariantEquations2D)

    _, _, _, rho_theta = u

    exner = equations.c_p * (rho_theta*equations.R/equations.p_0)^(equations.R/equations.c_v)
    return exner
end

end # @muladd
