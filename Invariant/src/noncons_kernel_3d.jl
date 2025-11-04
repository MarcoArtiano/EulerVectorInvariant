using Trixi
using Trixi: True, get_contravariant_vector, multiply_add_to_node_vars!, @threaded, get_surface_node_vars, get_normal_direction, get_node_coords, @turbo, PtrArray, StrideArray
using Trixi: StaticInt, indices
@muladd begin

	@inline function Trixi.flux_differencing_kernel!(du, u,
		element,
		mesh::Union{StructuredMesh{3}, P4estMesh{3},
			T8codeMesh{3}},
		nonconservative_terms::True, equations::CompressibleEulerVectorInvariantEquations3D,
		volume_flux, dg::DGSEM, cache, alpha = true)
		@unpack derivative_split = dg.basis
		@unpack contravariant_vectors = cache.elements
		symmetric_flux, nonconservative_flux = volume_flux
		# Apply the symmetric flux as usual
		Trixi.flux_differencing_kernel!(du, u, element, mesh, True(), equations, symmetric_flux,
			dg, cache, alpha)

		return nothing
	end

	# Inlined function for interface flux computation for flux + nonconservative terms
	@inline function Trixi.calc_interface_flux!(surface_flux_values,
		mesh::Union{P4estMesh{3}, T8codeMesh{3}},
		nonconservative_terms::True, equations::CompressibleEulerVectorInvariantEquations3D,
		surface_integral, dg::DG, cache,
		interface_index, normal_direction,
		primary_i_node_index, primary_j_node_index,
		primary_direction_index, primary_element_index,
		secondary_i_node_index, secondary_j_node_index,
		secondary_direction_index,
		secondary_element_index)
		@unpack u = cache.interfaces
		surface_flux, nonconservative_flux = surface_integral.surface_flux
		u_ll, u_rr = get_surface_node_vars(u, equations, dg, primary_i_node_index,
			primary_j_node_index, interface_index)

		flux_left, flux_right = surface_flux(u_ll, u_rr, normal_direction, equations)

		# Store the flux with nonconservative terms on the primary and secondary elements
		for v in eachvariable(equations)
			# Note the factor 0.5 necessary for the nonconservative fluxes based on
			# the interpretation of global SBP operators coupled discontinuously via
			# central fluxes/SATs
			surface_flux_values[v, primary_i_node_index, primary_j_node_index,
				primary_direction_index, primary_element_index] = flux_left[v]
			surface_flux_values[v, secondary_i_node_index, secondary_j_node_index,
				secondary_direction_index, secondary_element_index] = -flux_right[v]
		end

		return nothing
	end

	# inlined version of the boundary flux calculation along a physical interface
	@inline function Trixi.calc_boundary_flux!(surface_flux_values, t, boundary_condition,
		mesh::Union{P4estMesh{3}, T8codeMesh{3}},
		nonconservative_terms::True, equations::CompressibleEulerVectorInvariantEquations3D,
		surface_integral, dg::DG, cache, i_index, j_index,
		k_index, i_node_index, j_node_index,
		direction_index,
		element_index, boundary_index)
		@unpack boundaries = cache
		@unpack node_coordinates, contravariant_vectors = cache.elements
		@unpack surface_flux = surface_integral
		# Extract solution data from boundary container
		u_inner = get_node_vars(boundaries.u, equations, dg, i_node_index, j_node_index,
			boundary_index)

		# Outward-pointing normal direction (not normalized)
		normal_direction = get_normal_direction(direction_index, contravariant_vectors,
			i_index, j_index, k_index, element_index)

		# Coordinates at boundary node
		x = get_node_coords(node_coordinates, equations, dg,
			i_index, j_index, k_index, element_index)

		# Call pointwise numerical flux functions for the conservative and nonconservative part
		# in the normal direction on the boundary
		flux = boundary_condition(u_inner, normal_direction, x, t,
			surface_flux, equations)

		# Copy flux to element storage in the correct orientation
		for v in eachvariable(equations)
			# Note the factor 0.5 necessary for the nonconservative fluxes based on
			# the interpretation of global SBP operators coupled discontinuously via
			# central fluxes/SATs
			surface_flux_values[v, i_node_index, j_node_index,
				direction_index, element_index] = flux[v]
		end

		return nothing
	end

	@inline function flux_zero(u_ll, u_rr, normal_direction::AbstractVector, equations::CompressibleEulerVectorInvariantEquations3D)
		return zero(u_ll)
	end

	@inline function flux_invariant_adv_turbo(
		u_ll,
		u_rr,
		orientation_or_normal_direction,
		equations,
	)
		flux_invariant(u_ll, u_rr, orientation_or_normal_direction, equations)
	end
end

@inline function Trixi.flux_differencing_kernel!(_du::PtrArray, u_cons::PtrArray,
	element,
	mesh::Union{StructuredMesh{3}, P4estMesh{3}},
	nonconservative_terms::True,
	equations::CompressibleEulerVectorInvariantEquations3D,
	volume_flux::typeof(flux_invariant_turbo),
	dg::DGSEM, cache, alpha)
	@unpack derivative_split = dg.basis
	@unpack contravariant_vectors = cache.elements
	# Create a temporary array that will be used to store the RHS with permuted
	# indices `[i, j, k, v]` to allow using SIMD instructions.
	# `StrideArray`s with purely static dimensions do not allocate on the heap.
	du = StrideArray{eltype(u_cons)}(undef,
		(ntuple(_ -> StaticInt(nnodes(dg)), ndims(mesh))...,
			StaticInt(nvariables(equations))))

	# Convert conserved to primitive variables on the given `element`.
	u_prim = StrideArray{eltype(u_cons)}(undef,
		(ntuple(_ -> StaticInt(nnodes(dg)),
				ndims(mesh))...,
			StaticInt(nvariables(equations) + 1)))

	@turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
		rho = u_cons[1, i, j, k, element]
		v1 = u_cons[2, i, j, k, element]
		v2 = u_cons[3, i, j, k, element]
		v3 = u_cons[4, i, j, k, element]
		rho_theta = u_cons[5, i, j, k, element]
		phi = u_cons[6, i, j, k, element]
		theta = rho_theta/rho	
    		exner = (rho_theta * equations.R / equations.p_0)^(equations.R / equations.c_v)
		u_prim[i, j, k, 1] = rho
		u_prim[i, j, k, 2] = v1
		u_prim[i, j, k, 3] = v2
		u_prim[i, j, k, 4] = v3
		u_prim[i, j, k, 5] = theta
		u_prim[i, j, k, 6] = exner 
		u_prim[i, j, k, 7] = phi
	end

	# x direction
	# At first, we create new temporary arrays with permuted memory layout to
	# allow using SIMD instructions along the first dimension (which is contiguous
	# in memory).
	du_permuted = StrideArray{eltype(u_cons)}(undef,
		(StaticInt(nnodes(dg)^2),
			StaticInt(nnodes(dg)),
			StaticInt(nvariables(equations))))

	u_prim_permuted = StrideArray{eltype(u_cons)}(undef,
		(StaticInt(nnodes(dg)^2),
			StaticInt(nnodes(dg)),
			StaticInt(nvariables(equations) + 1)))

	@turbo for v in indices(u_prim, 4),
		k in eachnode(dg),
		j in eachnode(dg),
		i in eachnode(dg)

		jk = j + nnodes(dg) * (k - 1)
		u_prim_permuted[jk, i, v] = u_prim[i, j, k, v]
	end
	fill!(du_permuted, zero(eltype(du_permuted)))

	# We must also permute the contravariant vectors.
	contravariant_vectors_x = StrideArray{eltype(contravariant_vectors)}(undef,
		(StaticInt(nnodes(dg)^2),
			StaticInt(nnodes(dg)),
			StaticInt(ndims(mesh))))

	@turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
		jk = j + nnodes(dg) * (k - 1)
		contravariant_vectors_x[jk, i, 1] = contravariant_vectors[1, 1, i, j, k, element]
		contravariant_vectors_x[jk, i, 2] = contravariant_vectors[2, 1, i, j, k, element]
		contravariant_vectors_x[jk, i, 3] = contravariant_vectors[3, 1, i, j, k, element]
	end

	# Next, we basically inline the volume flux. To allow SIMD vectorization and
	# still use the symmetry of the volume flux and the derivative matrix, we
	# loop over the triangular part in an outer loop and use a plain inner loop.
	for i in eachnode(dg), ii in (i+1):nnodes(dg)
		@turbo for jk in Base.OneTo(nnodes(dg)^2)
			rho_ll = u_prim_permuted[jk, i, 1]
			v1_ll = u_prim_permuted[jk, i, 2]
			v2_ll = u_prim_permuted[jk, i, 3]
			v3_ll = u_prim_permuted[jk, i, 4]
			theta_ll = u_prim_permuted[jk, i, 5]
			exner_ll = u_prim_permuted[jk, i, 6]
			phi_ll = u_prim_permuted[jk, i, 7]

			rho_rr = u_prim_permuted[jk, ii, 1]
			v1_rr = u_prim_permuted[jk, ii, 2]
			v2_rr = u_prim_permuted[jk, ii, 3]
			v3_rr = u_prim_permuted[jk, ii, 4]
			theta_rr = u_prim_permuted[jk, ii, 5]
			exner_rr = u_prim_permuted[jk, ii, 6]
			phi_rr = u_prim_permuted[jk, ii, 7]

			normal_direction_1 = 0.5f0 * (contravariant_vectors_x[jk, i, 1] +
										contravariant_vectors_x[jk, ii, 1])
			normal_direction_2 = 0.5f0 * (contravariant_vectors_x[jk, i, 2] +
										contravariant_vectors_x[jk, ii, 2])
			normal_direction_3 = 0.5f0 * (contravariant_vectors_x[jk, i, 3] +
										contravariant_vectors_x[jk, ii, 3])

			v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 +
						 v3_ll * normal_direction_3
			v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 +
						 v3_rr * normal_direction_3


    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    kin_avg = 0.5f0 * (v1_rr * v1_rr + v2_rr * v2_rr + v3_rr * v3_rr + v1_ll * v1_ll + v2_ll * v2_ll + v3_ll * v3_ll)
    theta_avg = 0.5f0 * (theta_ll + theta_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)

    jump_v1 = v1_rr - v1_ll
    jump_v2 = v2_rr - v2_ll
    jump_v3 = v3_rr - v3_ll
	gravity = phi_rr - phi_ll
    f1 = rho_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = kin_avg * 0.5f0 * normal_direction_1 
    f3 = kin_avg * 0.5f0 * normal_direction_2 
    f4 = kin_avg * 0.5f0 * normal_direction_3 
    f5 = f1 * theta_avg
    theta_grad_exner =  equations.c_p * theta_avg * (exner_rr - exner_ll)    
    vorticity_x = v2_avg * (jump_v1 * normal_direction_2 - jump_v2 * normal_direction_1) + v3_avg * (jump_v1 * normal_direction_3 - jump_v3 * normal_direction_1)
    vorticity_y = v1_avg * (jump_v2 * normal_direction_1 - jump_v1 * normal_direction_2) + v3_avg * (jump_v2 * normal_direction_3 - jump_v3 * normal_direction_2)
    vorticity_z = v1_avg * (jump_v3 * normal_direction_1 - jump_v1 * normal_direction_3) + v2_avg * (jump_v3 * normal_direction_2 - jump_v2 * normal_direction_3) 
 
			g2 = vorticity_x + (theta_grad_exner + gravity)* normal_direction_1
			g3 = vorticity_y + (theta_grad_exner + gravity)* normal_direction_2
			g4 = vorticity_z + (theta_grad_exner + gravity) * normal_direction_3

			# Add scaled fluxes to RHS
			factor_i = alpha * derivative_split[i, ii]
			du_permuted[jk, i, 1] += factor_i * f1
			du_permuted[jk, i, 2] += factor_i * (f2 + 0.5 * g2)
			du_permuted[jk, i, 3] += factor_i * (f3 + 0.5 * g3)
			du_permuted[jk, i, 4] += factor_i * (f4 + 0.5 * g4)
			du_permuted[jk, i, 5] += factor_i * f5

			factor_ii = alpha * derivative_split[ii, i]
			du_permuted[jk, ii, 1] += factor_ii * f1
			du_permuted[jk, ii, 2] += factor_ii * (f2 - 0.5 * g2)
			du_permuted[jk, ii, 3] += factor_ii * (f3 - 0.5 * g3)
			du_permuted[jk, ii, 4] += factor_ii * (f4 - 0.5 * g4)
			du_permuted[jk, ii, 5] += factor_ii * f5
		end
	end

	@turbo for v in eachvariable(equations),
		k in eachnode(dg),
		j in eachnode(dg),
		i in eachnode(dg)

		jk = j + nnodes(dg) * (k - 1)
		du[i, j, k, v] = du_permuted[jk, i, v]
	end

	# y direction
	# We must also permute the contravariant vectors.
	contravariant_vectors_y = StrideArray{eltype(contravariant_vectors)}(undef,
		(StaticInt(nnodes(dg)),
			StaticInt(nnodes(dg)),
			StaticInt(nnodes(dg)),
			StaticInt(ndims(mesh))))

	@turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
		contravariant_vectors_y[i, j, k, 1] = contravariant_vectors[1, 2, i, j, k, element]
		contravariant_vectors_y[i, j, k, 2] = contravariant_vectors[2, 2, i, j, k, element]
		contravariant_vectors_y[i, j, k, 3] = contravariant_vectors[3, 2, i, j, k, element]
	end

	# A possible permutation of array dimensions with improved opportunities for
	# SIMD vectorization appeared to be slower than the direct version used here
	# in preliminary numerical experiments on an AVX2 system.
	for j in eachnode(dg), jj in (j+1):nnodes(dg)
		@turbo for k in eachnode(dg), i in eachnode(dg)
			rho_ll = u_prim[i, j, k, 1]
			v1_ll = u_prim[i, j, k, 2]
			v2_ll = u_prim[i, j, k, 3]
			v3_ll = u_prim[i, j, k, 4]
			theta_ll = u_prim[i, j, k, 5]
			exner_ll = u_prim[i, j, k, 6]
			phi_ll = u_prim[i, j, k, 7]

			rho_rr = u_prim[i, jj, k, 1]
			v1_rr = u_prim[i, jj, k, 2]
			v2_rr = u_prim[i, jj, k, 3]
			v3_rr = u_prim[i, jj, k, 4]
			theta_rr = u_prim[i, jj, k, 5]
			exner_rr = u_prim[i, jj, k, 6]
			phi_rr = u_prim[i, jj, k, 7]

			normal_direction_1 = 0.5f0 * (contravariant_vectors_y[i, j, k, 1] +
										contravariant_vectors_y[i, jj, k, 1])
			normal_direction_2 = 0.5f0 * (contravariant_vectors_y[i, j, k, 2] +
										contravariant_vectors_y[i, jj, k, 2])
			normal_direction_3 = 0.5f0 * (contravariant_vectors_y[i, j, k, 3] +
										contravariant_vectors_y[i, jj, k, 3])

			v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 +
						 v3_ll * normal_direction_3
			v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 +
						 v3_rr * normal_direction_3
			
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    kin_avg = 0.5f0 * (v1_rr * v1_rr + v2_rr * v2_rr + v3_rr * v3_rr + v1_ll * v1_ll + v2_ll * v2_ll + v3_ll * v3_ll)
    theta_avg = 0.5f0 * (theta_ll + theta_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)

    jump_v1 = v1_rr - v1_ll
    jump_v2 = v2_rr - v2_ll
    jump_v3 = v3_rr - v3_ll
    gravity = phi_rr - phi_ll
    f1 = rho_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = kin_avg * 0.5f0 * normal_direction_1 
    f3 = kin_avg * 0.5f0 * normal_direction_2 
    f4 = kin_avg * 0.5f0 * normal_direction_3 
    f5 = f1 * theta_avg
    theta_grad_exner =  equations.c_p * theta_avg * (exner_rr - exner_ll)    
    vorticity_x = v2_avg * (jump_v1 * normal_direction_2 - jump_v2 * normal_direction_1) + v3_avg * (jump_v1 * normal_direction_3 - jump_v3 * normal_direction_1)
    vorticity_y = v1_avg * (jump_v2 * normal_direction_1 - jump_v1 * normal_direction_2) + v3_avg * (jump_v2 * normal_direction_3 - jump_v3 * normal_direction_2)
    vorticity_z = v1_avg * (jump_v3 * normal_direction_1 - jump_v1 * normal_direction_3) + v2_avg * (jump_v3 * normal_direction_2 - jump_v2 * normal_direction_3) 
 
			g2 = vorticity_x + (theta_grad_exner + gravity)* normal_direction_1
			g3 = vorticity_y + (theta_grad_exner + gravity)* normal_direction_2
			g4 = vorticity_z + (theta_grad_exner + gravity) * normal_direction_3

			# Add scaled fluxes to RHS
			factor_j = alpha * derivative_split[j, jj]
			du[i, j, k, 1] += factor_j * f1
			du[i, j, k, 2] += factor_j * (f2 + 0.5 * g2)
			du[i, j, k, 3] += factor_j * (f3 + 0.5 * g3)
			du[i, j, k, 4] += factor_j * (f4 + 0.5 * g4)
			du[i, j, k, 5] += factor_j * f5

			factor_jj = alpha * derivative_split[jj, j]
			du[i, jj, k, 1] += factor_jj * f1
			du[i, jj, k, 2] += factor_jj * (f2 - 0.5 * g2)
			du[i, jj, k, 3] += factor_jj * (f3 - 0.5 * g3)
			du[i, jj, k, 4] += factor_jj * (f4 - 0.5 * g4)
			du[i, jj, k, 5] += factor_jj * f5
		end
	end

	# z direction
	# The memory layout is already optimal for SIMD vectorization in this loop.
	# We just squeeze the first two dimensions to make the code slightly faster.
	GC.@preserve u_prim begin
		u_prim_reshaped = PtrArray(pointer(u_prim),
			(StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
				StaticInt(nvariables(equations) + 1)))

		du_reshaped = PtrArray(pointer(du),
			(StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
				StaticInt(nvariables(equations))))

		# We must also permute the contravariant vectors.
		contravariant_vectors_z = StrideArray{eltype(contravariant_vectors)}(undef,
			(StaticInt(nnodes(dg)^2),
				StaticInt(nnodes(dg)),
				StaticInt(ndims(mesh))))

		@turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
			ij = i + nnodes(dg) * (j - 1)
			contravariant_vectors_z[ij, k, 1] = contravariant_vectors[1, 3, i, j, k,
				element]
			contravariant_vectors_z[ij, k, 2] = contravariant_vectors[2, 3, i, j, k,
				element]
			contravariant_vectors_z[ij, k, 3] = contravariant_vectors[3, 3, i, j, k,
				element]
		end

		for k in eachnode(dg), kk in (k+1):nnodes(dg)
			@turbo for ij in Base.OneTo(nnodes(dg)^2)
				rho_ll = u_prim_reshaped[ij, k, 1]
				v1_ll = u_prim_reshaped[ij, k, 2]
				v2_ll = u_prim_reshaped[ij, k, 3]
				v3_ll = u_prim_reshaped[ij, k, 4]
				theta_ll = u_prim_reshaped[ij, k, 5]
				exner_ll = u_prim_reshaped[ij, k, 6]
				phi_ll = u_prim_reshaped[ij, k, 7]


				rho_rr = u_prim_reshaped[ij, kk, 1]
				v1_rr = u_prim_reshaped[ij, kk, 2]
				v2_rr = u_prim_reshaped[ij, kk, 3]
				v3_rr = u_prim_reshaped[ij, kk, 4]
				theta_rr = u_prim_reshaped[ij, kk, 5]
				exner_rr = u_prim_reshaped[ij, kk, 6]
				phi_rr = u_prim_reshaped[ij, kk, 7]


				normal_direction_1 = 0.5f0 * (contravariant_vectors_z[ij, k, 1] +
											contravariant_vectors_z[ij, kk, 1])
				normal_direction_2 = 0.5f0 * (contravariant_vectors_z[ij, k, 2] +
											contravariant_vectors_z[ij, kk, 2])
				normal_direction_3 = 0.5f0 * (contravariant_vectors_z[ij, k, 3] +
											contravariant_vectors_z[ij, kk, 3])

				v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 +
							 v3_ll * normal_direction_3
				v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 +
							 v3_rr * normal_direction_3

    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    kin_avg = 0.5f0 * (v1_rr * v1_rr + v2_rr * v2_rr + v3_rr * v3_rr + v1_ll * v1_ll + v2_ll * v2_ll + v3_ll * v3_ll)
    theta_avg = 0.5f0 * (theta_ll + theta_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)

    jump_v1 = v1_rr - v1_ll
    jump_v2 = v2_rr - v2_ll
    jump_v3 = v3_rr - v3_ll
    gravity = phi_rr - phi_ll
    f1 = rho_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = kin_avg * 0.5f0 * normal_direction_1 
    f3 = kin_avg * 0.5f0 * normal_direction_2 
    f4 = kin_avg * 0.5f0 * normal_direction_3 
    f5 = f1 * theta_avg
    theta_grad_exner =  equations.c_p * theta_avg * (exner_rr - exner_ll)    
    vorticity_x = v2_avg * (jump_v1 * normal_direction_2- jump_v2 * normal_direction_1) + v3_avg * (jump_v1 * normal_direction_3 - jump_v3 * normal_direction_1)
    vorticity_y = v1_avg * (jump_v2 * normal_direction_1 - jump_v1 * normal_direction_2) + v3_avg * (jump_v2 * normal_direction_3 - jump_v3 * normal_direction_2)
    vorticity_z = v1_avg * (jump_v3 * normal_direction_1 - jump_v1 * normal_direction_3) + v2_avg * (jump_v3 * normal_direction_2 - jump_v2 * normal_direction_3) 
 
			g2 = vorticity_x + (theta_grad_exner + gravity)* normal_direction_1
			g3 = vorticity_y + (theta_grad_exner + gravity)* normal_direction_2
			g4 = vorticity_z + (theta_grad_exner + gravity) * normal_direction_3
				# Add scaled fluxes to RHS
				factor_k = alpha * derivative_split[k, kk]
				du_reshaped[ij, k, 1] += factor_k * f1
				du_reshaped[ij, k, 2] += factor_k * (f2 + 0.5 * g2)
				du_reshaped[ij, k, 3] += factor_k * (f3 + 0.5 * g3)
				du_reshaped[ij, k, 4] += factor_k * (f4 + 0.5 * g4)
				du_reshaped[ij, k, 5] += factor_k * f5

				factor_kk = alpha * derivative_split[kk, k]
				du_reshaped[ij, kk, 1] += factor_kk * f1
				du_reshaped[ij, kk, 2] += factor_kk * (f2 - 0.5 * g2)
				du_reshaped[ij, kk, 3] += factor_kk * (f3 - 0.5 * g3)
				du_reshaped[ij, kk, 4] += factor_kk * (f4 - 0.5 * g4)
				du_reshaped[ij, kk, 5] += factor_kk * f5
			end
		end
	end # GC.@preserve u_prim begin

	# Finally, we add the temporary RHS computed here to the global RHS in the
	# given `element`.
	@turbo for v in eachvariable(equations),
		k in eachnode(dg),
		j in eachnode(dg),
		i in eachnode(dg)

		_du[v, i, j, k, element] += du[i, j, k, v]
	end
end

@inline function Trixi.flux_differencing_kernel!(_du::PtrArray, u_cons::PtrArray,
	element,
	mesh::Union{StructuredMesh{3}, P4estMesh{3}},
	nonconservative_terms::True,
	equations::CompressibleEulerVectorInvariantEquations3D,
	volume_flux::typeof(flux_invariant_adv_turbo),
	dg::DGSEM, cache, alpha)
	@unpack derivative_split = dg.basis
	@unpack contravariant_vectors = cache.elements
	# Create a temporary array that will be used to store the RHS with permuted
	# indices `[i, j, k, v]` to allow using SIMD instructions.
	# `StrideArray`s with purely static dimensions do not allocate on the heap.
	du = StrideArray{eltype(u_cons)}(undef,
		(ntuple(_ -> StaticInt(nnodes(dg)), ndims(mesh))...,
			StaticInt(nvariables(equations))))

	# Convert conserved to primitive variables on the given `element`.
	u_prim = StrideArray{eltype(u_cons)}(undef,
		(ntuple(_ -> StaticInt(nnodes(dg)),
				ndims(mesh))...,
			StaticInt(nvariables(equations) + 1)))

	@turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
		rho = u_cons[1, i, j, k, element]
		v1 = u_cons[2, i, j, k, element]
		v2 = u_cons[3, i, j, k, element]
		v3 = u_cons[4, i, j, k, element]
		rho_theta = u_cons[5, i, j, k, element]
		phi = u_cons[6, i, j, k, element]
		theta = rho_theta/rho	
    		exner = (rho_theta * equations.R / equations.p_0)^(equations.R / equations.c_v)
		u_prim[i, j, k, 1] = rho
		u_prim[i, j, k, 2] = v1
		u_prim[i, j, k, 3] = v2
		u_prim[i, j, k, 4] = v3
		u_prim[i, j, k, 5] = theta
		u_prim[i, j, k, 6] = exner 
		u_prim[i, j, k, 7] = phi
	end

	# x direction
	# At first, we create new temporary arrays with permuted memory layout to
	# allow using SIMD instructions along the first dimension (which is contiguous
	# in memory).
	du_permuted = StrideArray{eltype(u_cons)}(undef,
		(StaticInt(nnodes(dg)^2),
			StaticInt(nnodes(dg)),
			StaticInt(nvariables(equations))))

	u_prim_permuted = StrideArray{eltype(u_cons)}(undef,
		(StaticInt(nnodes(dg)^2),
			StaticInt(nnodes(dg)),
			StaticInt(nvariables(equations) + 1)))

	@turbo for v in indices(u_prim, 4),
		k in eachnode(dg),
		j in eachnode(dg),
		i in eachnode(dg)

		jk = j + nnodes(dg) * (k - 1)
		u_prim_permuted[jk, i, v] = u_prim[i, j, k, v]
	end
	fill!(du_permuted, zero(eltype(du_permuted)))

	# We must also permute the contravariant vectors.
	contravariant_vectors_x = StrideArray{eltype(contravariant_vectors)}(undef,
		(StaticInt(nnodes(dg)^2),
			StaticInt(nnodes(dg)),
			StaticInt(ndims(mesh))))

	@turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
		jk = j + nnodes(dg) * (k - 1)
		contravariant_vectors_x[jk, i, 1] = contravariant_vectors[1, 1, i, j, k, element]
		contravariant_vectors_x[jk, i, 2] = contravariant_vectors[2, 1, i, j, k, element]
		contravariant_vectors_x[jk, i, 3] = contravariant_vectors[3, 1, i, j, k, element]
	end

	# Next, we basically inline the volume flux. To allow SIMD vectorization and
	# still use the symmetry of the volume flux and the derivative matrix, we
	# loop over the triangular part in an outer loop and use a plain inner loop.
	for i in eachnode(dg), ii in (i+1):nnodes(dg)
		@turbo for jk in Base.OneTo(nnodes(dg)^2)
			rho_ll = u_prim_permuted[jk, i, 1]
			v1_ll = u_prim_permuted[jk, i, 2]
			v2_ll = u_prim_permuted[jk, i, 3]
			v3_ll = u_prim_permuted[jk, i, 4]
			theta_ll = u_prim_permuted[jk, i, 5]
			exner_ll = u_prim_permuted[jk, i, 6]
			phi_ll = u_prim_permuted[jk, i, 7]

			rho_rr = u_prim_permuted[jk, ii, 1]
			v1_rr = u_prim_permuted[jk, ii, 2]
			v2_rr = u_prim_permuted[jk, ii, 3]
			v3_rr = u_prim_permuted[jk, ii, 4]
			theta_rr = u_prim_permuted[jk, ii, 5]
			exner_rr = u_prim_permuted[jk, ii, 6]
			phi_rr = u_prim_permuted[jk, ii, 7]

			normal_direction_1 = 0.5f0 * (contravariant_vectors_x[jk, i, 1] +
										contravariant_vectors_x[jk, ii, 1])
			normal_direction_2 = 0.5f0 * (contravariant_vectors_x[jk, i, 2] +
										contravariant_vectors_x[jk, ii, 2])
			normal_direction_3 = 0.5f0 * (contravariant_vectors_x[jk, i, 3] +
										contravariant_vectors_x[jk, ii, 3])

			v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 +
						 v3_ll * normal_direction_3
			v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 +
						 v3_rr * normal_direction_3


    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    theta_avg = 0.5f0 * (theta_ll + theta_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)

    jump_v1 = v1_rr - v1_ll
    jump_v2 = v2_rr - v2_ll
    jump_v3 = v3_rr - v3_ll
    gravity = phi_rr - phi_ll
    f1 = rho_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = 0
    f3 = 0 
    f4 = 0 
    f5 = f1 * theta_avg
    theta_grad_exner =  equations.c_p * theta_avg * (exner_rr - exner_ll)    
 
    advection = v1_avg * normal_direction_1 + v2_avg * normal_direction_2 + v3_avg * normal_direction_3 
    advection_x = advection * jump_v1
    advection_y = advection * jump_v2
    advection_z = advection * jump_v3
			g2 = advection_x + (theta_grad_exner + gravity)* normal_direction_1
			g3 = advection_y + (theta_grad_exner + gravity)* normal_direction_2
			g4 = advection_z + (theta_grad_exner + gravity) * normal_direction_3

			# Add scaled fluxes to RHS
			factor_i = alpha * derivative_split[i, ii]
			du_permuted[jk, i, 1] += factor_i * f1
			du_permuted[jk, i, 2] += factor_i * (f2 + 0.5 * g2)
			du_permuted[jk, i, 3] += factor_i * (f3 + 0.5 * g3)
			du_permuted[jk, i, 4] += factor_i * (f4 + 0.5 * g4)
			du_permuted[jk, i, 5] += factor_i * f5

			factor_ii = alpha * derivative_split[ii, i]
			du_permuted[jk, ii, 1] += factor_ii * f1
			du_permuted[jk, ii, 2] += factor_ii * (f2 - 0.5 * g2)
			du_permuted[jk, ii, 3] += factor_ii * (f3 - 0.5 * g3)
			du_permuted[jk, ii, 4] += factor_ii * (f4 - 0.5 * g4)
			du_permuted[jk, ii, 5] += factor_ii * f5
		end
	end

	@turbo for v in eachvariable(equations),
		k in eachnode(dg),
		j in eachnode(dg),
		i in eachnode(dg)

		jk = j + nnodes(dg) * (k - 1)
		du[i, j, k, v] = du_permuted[jk, i, v]
	end

	# y direction
	# We must also permute the contravariant vectors.
	contravariant_vectors_y = StrideArray{eltype(contravariant_vectors)}(undef,
		(StaticInt(nnodes(dg)),
			StaticInt(nnodes(dg)),
			StaticInt(nnodes(dg)),
			StaticInt(ndims(mesh))))

	@turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
		contravariant_vectors_y[i, j, k, 1] = contravariant_vectors[1, 2, i, j, k, element]
		contravariant_vectors_y[i, j, k, 2] = contravariant_vectors[2, 2, i, j, k, element]
		contravariant_vectors_y[i, j, k, 3] = contravariant_vectors[3, 2, i, j, k, element]
	end

	# A possible permutation of array dimensions with improved opportunities for
	# SIMD vectorization appeared to be slower than the direct version used here
	# in preliminary numerical experiments on an AVX2 system.
	for j in eachnode(dg), jj in (j+1):nnodes(dg)
		@turbo for k in eachnode(dg), i in eachnode(dg)
			rho_ll = u_prim[i, j, k, 1]
			v1_ll = u_prim[i, j, k, 2]
			v2_ll = u_prim[i, j, k, 3]
			v3_ll = u_prim[i, j, k, 4]
			theta_ll = u_prim[i, j, k, 5]
			exner_ll = u_prim[i, j, k, 6]
			phi_ll = u_prim[i, j, k, 7]

			rho_rr = u_prim[i, jj, k, 1]
			v1_rr = u_prim[i, jj, k, 2]
			v2_rr = u_prim[i, jj, k, 3]
			v3_rr = u_prim[i, jj, k, 4]
			theta_rr = u_prim[i, jj, k, 5]
			exner_rr = u_prim[i, jj, k, 6]
			phi_rr = u_prim[i, jj, k, 7]

			normal_direction_1 = 0.5f0 * (contravariant_vectors_y[i, j, k, 1] +
										contravariant_vectors_y[i, jj, k, 1])
			normal_direction_2 = 0.5f0 * (contravariant_vectors_y[i, j, k, 2] +
										contravariant_vectors_y[i, jj, k, 2])
			normal_direction_3 = 0.5f0 * (contravariant_vectors_y[i, j, k, 3] +
										contravariant_vectors_y[i, jj, k, 3])

			v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 +
						 v3_ll * normal_direction_3
			v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 +
						 v3_rr * normal_direction_3
			
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    kin_avg = 0.5f0 * (v1_rr * v1_rr + v2_rr * v2_rr + v3_rr * v3_rr + v1_ll * v1_ll + v2_ll * v2_ll + v3_ll * v3_ll)
    theta_avg = 0.5f0 * (theta_ll + theta_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)

    jump_v1 = v1_rr - v1_ll
    jump_v2 = v2_rr - v2_ll
    jump_v3 = v3_rr - v3_ll
    gravity = phi_rr - phi_ll
    f1 = rho_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = 0
    f3 = 0 
    f4 = 0 
    f5 = f1 * theta_avg
    theta_grad_exner =  equations.c_p * theta_avg * (exner_rr - exner_ll)    
 
    advection = v1_avg * normal_direction_1 + v2_avg * normal_direction_2 + v3_avg * normal_direction_3 
    advection_x = advection * jump_v1
    advection_y = advection * jump_v2
    advection_z = advection * jump_v3
			g2 = advection_x + (theta_grad_exner + gravity)* normal_direction_1
			g3 = advection_y + (theta_grad_exner + gravity)* normal_direction_2
			g4 = advection_z + (theta_grad_exner + gravity) * normal_direction_3

			# Add scaled fluxes to RHS
			factor_j = alpha * derivative_split[j, jj]
			du[i, j, k, 1] += factor_j * f1
			du[i, j, k, 2] += factor_j * (f2 + 0.5 * g2)
			du[i, j, k, 3] += factor_j * (f3 + 0.5 * g3)
			du[i, j, k, 4] += factor_j * (f4 + 0.5 * g4)
			du[i, j, k, 5] += factor_j * f5

			factor_jj = alpha * derivative_split[jj, j]
			du[i, jj, k, 1] += factor_jj * f1
			du[i, jj, k, 2] += factor_jj * (f2 - 0.5 * g2)
			du[i, jj, k, 3] += factor_jj * (f3 - 0.5 * g3)
			du[i, jj, k, 4] += factor_jj * (f4 - 0.5 * g4)
			du[i, jj, k, 5] += factor_jj * f5
		end
	end

	# z direction
	# The memory layout is already optimal for SIMD vectorization in this loop.
	# We just squeeze the first two dimensions to make the code slightly faster.
	GC.@preserve u_prim begin
		u_prim_reshaped = PtrArray(pointer(u_prim),
			(StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
				StaticInt(nvariables(equations) + 1)))

		du_reshaped = PtrArray(pointer(du),
			(StaticInt(nnodes(dg)^2), StaticInt(nnodes(dg)),
				StaticInt(nvariables(equations))))

		# We must also permute the contravariant vectors.
		contravariant_vectors_z = StrideArray{eltype(contravariant_vectors)}(undef,
			(StaticInt(nnodes(dg)^2),
				StaticInt(nnodes(dg)),
				StaticInt(ndims(mesh))))

		@turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
			ij = i + nnodes(dg) * (j - 1)
			contravariant_vectors_z[ij, k, 1] = contravariant_vectors[1, 3, i, j, k,
				element]
			contravariant_vectors_z[ij, k, 2] = contravariant_vectors[2, 3, i, j, k,
				element]
			contravariant_vectors_z[ij, k, 3] = contravariant_vectors[3, 3, i, j, k,
				element]
		end

		for k in eachnode(dg), kk in (k+1):nnodes(dg)
			@turbo for ij in Base.OneTo(nnodes(dg)^2)
				rho_ll = u_prim_reshaped[ij, k, 1]
				v1_ll = u_prim_reshaped[ij, k, 2]
				v2_ll = u_prim_reshaped[ij, k, 3]
				v3_ll = u_prim_reshaped[ij, k, 4]
				theta_ll = u_prim_reshaped[ij, k, 5]
				exner_ll = u_prim_reshaped[ij, k, 6]
				phi_ll = u_prim_reshaped[ij, k, 7]


				rho_rr = u_prim_reshaped[ij, kk, 1]
				v1_rr = u_prim_reshaped[ij, kk, 2]
				v2_rr = u_prim_reshaped[ij, kk, 3]
				v3_rr = u_prim_reshaped[ij, kk, 4]
				theta_rr = u_prim_reshaped[ij, kk, 5]
				exner_rr = u_prim_reshaped[ij, kk, 6]
				phi_rr = u_prim_reshaped[ij, kk, 7]


				normal_direction_1 = 0.5f0 * (contravariant_vectors_z[ij, k, 1] +
											contravariant_vectors_z[ij, kk, 1])
				normal_direction_2 = 0.5f0 * (contravariant_vectors_z[ij, k, 2] +
											contravariant_vectors_z[ij, kk, 2])
				normal_direction_3 = 0.5f0 * (contravariant_vectors_z[ij, k, 3] +
											contravariant_vectors_z[ij, kk, 3])

				v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2 +
							 v3_ll * normal_direction_3
				v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2 +
							 v3_rr * normal_direction_3

    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    kin_avg = 0.5f0 * (v1_rr * v1_rr + v2_rr * v2_rr + v3_rr * v3_rr + v1_ll * v1_ll + v2_ll * v2_ll + v3_ll * v3_ll)
    theta_avg = 0.5f0 * (theta_ll + theta_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)

    jump_v1 = v1_rr - v1_ll
    jump_v2 = v2_rr - v2_ll
    jump_v3 = v3_rr - v3_ll
    gravity = phi_rr - phi_ll

    f1 = rho_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = 0
    f3 = 0 
    f4 = 0 
    f5 = f1 * theta_avg
    theta_grad_exner =  equations.c_p * theta_avg * (exner_rr - exner_ll)    
 
    advection = v1_avg * normal_direction_1 + v2_avg * normal_direction_2 + v3_avg * normal_direction_3 
    advection_x = advection * jump_v1
    advection_y = advection * jump_v2
    advection_z = advection * jump_v3
			g2 = advection_x + (theta_grad_exner + gravity)* normal_direction_1
			g3 = advection_y + (theta_grad_exner + gravity)* normal_direction_2
			g4 = advection_z + (theta_grad_exner + gravity) * normal_direction_3
				# Add scaled fluxes to RHS
				factor_k = alpha * derivative_split[k, kk]
				du_reshaped[ij, k, 1] += factor_k * f1
				du_reshaped[ij, k, 2] += factor_k * (f2 + 0.5 * g2)
				du_reshaped[ij, k, 3] += factor_k * (f3 + 0.5 * g3)
				du_reshaped[ij, k, 4] += factor_k * (f4 + 0.5 * g4)
				du_reshaped[ij, k, 5] += factor_k * f5

				factor_kk = alpha * derivative_split[kk, k]
				du_reshaped[ij, kk, 1] += factor_kk * f1
				du_reshaped[ij, kk, 2] += factor_kk * (f2 - 0.5 * g2)
				du_reshaped[ij, kk, 3] += factor_kk * (f3 - 0.5 * g3)
				du_reshaped[ij, kk, 4] += factor_kk * (f4 - 0.5 * g4)
				du_reshaped[ij, kk, 5] += factor_kk * f5
			end
		end
	end # GC.@preserve u_prim begin

	# Finally, we add the temporary RHS computed here to the global RHS in the
	# given `element`.
	@turbo for v in eachvariable(equations),
		k in eachnode(dg),
		j in eachnode(dg),
		i in eachnode(dg)

		_du[v, i, j, k, element] += du[i, j, k, v]
	end
end
