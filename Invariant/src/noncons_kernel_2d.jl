using Trixi
using Trixi: True, get_contravariant_vector, multiply_add_to_node_vars!, @threaded, get_surface_node_vars, get_normal_direction, get_node_coords, @turbo, PtrArray, StrideArray
using Trixi: StaticInt, indices
@muladd begin

	@inline function Trixi.flux_differencing_kernel!(du, u,
		element,
		mesh::Union{StructuredMesh{2},
			StructuredMeshView{2},
			UnstructuredMesh2D, P4estMesh{2},
			T8codeMesh{2}},
		nonconservative_terms::True, equations::CompressibleEulerVectorInvariantEquations2D,
		volume_flux, dg::DGSEM, cache, alpha = true)
		@unpack derivative_split = dg.basis
		@unpack contravariant_vectors = cache.elements
		symmetric_flux, nonconservative_flux = volume_flux

		Trixi.flux_differencing_kernel!(du, u, element, mesh, True(), equations, symmetric_flux,
			dg, cache, alpha)

		return nothing
	end

	# Inlined version of the interface flux computation for equations with conservative and nonconservative terms
	@inline function Trixi.calc_interface_flux!(surface_flux_values,
		mesh::Union{P4estMesh{2}, T8codeMesh{2}},
		nonconservative_terms::True, equations::CompressibleEulerVectorInvariantEquations2D,
		surface_integral, dg::DG, cache,
		interface_index, normal_direction,
		primary_node_index, primary_direction_index,
		primary_element_index,
		secondary_node_index, secondary_direction_index,
		secondary_element_index)
		@unpack u = cache.interfaces
		surface_flux, nonconservative_flux = surface_integral.surface_flux

		u_ll, u_rr = get_surface_node_vars(u, equations, dg, primary_node_index,
			interface_index)

		flux_left, flux_right = surface_flux(u_ll, u_rr, normal_direction, equations)

		# Store the flux with nonconservative terms on the primary and secondary elements
		for v in eachvariable(equations)
			# Note the factor 0.5 necessary for the nonconservative fluxes based on
			# the interpretation of global SBP operators coupled discontinuously via
			# central fluxes/SATs
			surface_flux_values[v, primary_node_index, primary_direction_index, primary_element_index] = flux_left[v]
			surface_flux_values[v, secondary_node_index, secondary_direction_index, secondary_element_index] = -flux_right[v]
		end

		return nothing
	end

	@inline function Trixi.calc_boundary_flux!(surface_flux_values, t, boundary_condition,
		mesh::Union{P4estMesh{2}, T8codeMesh{2}},
		nonconservative_terms::True, equations::CompressibleEulerVectorInvariantEquations2D,
		surface_integral, dg::DG, cache,
		i_index, j_index,
		node_index, direction_index, element_index,
		boundary_index)
		@unpack boundaries = cache
		@unpack node_coordinates, contravariant_vectors = cache.elements

		# Extract solution data from boundary container
		u_inner = get_node_vars(boundaries.u, equations, dg, node_index, boundary_index)

		# Outward-pointing normal direction (not normalized)
		normal_direction = get_normal_direction(direction_index, contravariant_vectors,
			i_index, j_index, element_index)

		# Coordinates at boundary node
		x = get_node_coords(node_coordinates, equations, dg,
			i_index, j_index, element_index)

		# Call pointwise numerical flux functions for the conservative and nonconservative part
		# in the normal direction on the boundary
		flux = boundary_condition(u_inner, normal_direction, x, t,
			surface_integral.surface_flux, equations)

		# Copy flux to element storage in the correct orientation
		for v in eachvariable(equations)
			# Note the factor 0.5 necessary for the nonconservative fluxes based on
			# the interpretation of global SBP operators coupled discontinuously via
			# central fluxes/SATs
			surface_flux_values[v, node_index, direction_index, element_index] = flux[v]
		end

		return nothing
	end

	@inline function flux_invariant_turbo(u_ll, u_rr, orientation_or_normal_direction,
		equations)
		flux_invariant(u_ll, u_rr, orientation_or_normal_direction, equations)
	end

	@inline function flux_zero(u_ll, u_rr, normal_direction::AbstractVector, equations)
	return zero(u_ll)
	end

end

@inline function Trixi.flux_differencing_kernel!(_du::PtrArray, u_cons::PtrArray,
	element,
	mesh::Union{StructuredMesh{2},
		UnstructuredMesh2D, P4estMesh{2}},
	nonconservative_terms::True,
	equations::CompressibleEulerVectorInvariantEquations2D,
	volume_flux::typeof(flux_invariant_turbo),
	dg::DGSEM, cache, alpha)
	@unpack derivative_split = dg.basis
	@unpack contravariant_vectors = cache.elements
	du = StrideArray{eltype(u_cons)}(undef,
		(ntuple(_ -> StaticInt(nnodes(dg)), ndims(mesh))...,
			StaticInt(nvariables(equations))))

	u_prim = StrideArray{eltype(u_cons)}(undef,
		(ntuple(_ -> StaticInt(nnodes(dg)),
				ndims(mesh))...,
			StaticInt(nvariables(equations) + 1))) # We also compute "+ 2" logs

	@turbo for j in eachnode(dg), i in eachnode(dg)
		rho = u_cons[1, i, j, element]
		rho_v1 = u_cons[2, i, j, element]
		rho_v2 = u_cons[3, i, j, element]
		rho_theta = u_cons[4, i, j, element]
		phi = u_cons[5, i, j, element]
		v1 = rho_v1 / rho
		v2 = rho_v2 / rho
		theta = rho_theta/rho
		exner = (rho_theta * equations.R / equations.p_0)^(equations.R / equations.c_v)
		u_prim[i, j, 1] = rho
		u_prim[i, j, 2] = v1
		u_prim[i, j, 3] = v2
		u_prim[i, j, 4] = theta
		u_prim[i, j, 5] = exner
		u_prim[i, j, 6] = phi
	end

	du_permuted = StrideArray{eltype(u_cons)}(undef,
		(StaticInt(nnodes(dg)), StaticInt(nnodes(dg)),
			StaticInt(nvariables(equations))))

	u_prim_permuted = StrideArray{eltype(u_cons)}(undef,
		(StaticInt(nnodes(dg)),
			StaticInt(nnodes(dg)),
			StaticInt(nvariables(equations) + 1)))

	@turbo for v in indices(u_prim, 3), # v in eachvariable(equations) misses +2 logs
		j in eachnode(dg),
		i in eachnode(dg)

		u_prim_permuted[j, i, v] = u_prim[i, j, v]
	end
	fill!(du_permuted, zero(eltype(du_permuted)))

	contravariant_vectors_x = StrideArray{eltype(contravariant_vectors)}(undef,
		(StaticInt(nnodes(dg)),
			StaticInt(nnodes(dg)),
			StaticInt(ndims(mesh))))

	@turbo for j in eachnode(dg), i in eachnode(dg)
		contravariant_vectors_x[j, i, 1] = contravariant_vectors[1, 1, i, j, element]
		contravariant_vectors_x[j, i, 2] = contravariant_vectors[2, 1, i, j, element]
	end

	# Next, we basically inline the volume flux. To allow SIMD vectorization and
	# still use the symmetry of the volume flux and the derivative matrix, we
	# loop over the triangular part in an outer loop and use a plain inner loop.
	for i in eachnode(dg), ii in (i+1):nnodes(dg)
		@turbo for j in eachnode(dg)
			rho_ll = u_prim_permuted[j, i, 1]
			v1_ll = u_prim_permuted[j, i, 2]
			v2_ll = u_prim_permuted[j, i, 3]
			theta_ll = u_prim_permuted[j, i, 4]
			exner_ll = u_prim_permuted[j, i, 5]
			phi_ll = u_prim_permuted[j, i, 6]

			rho_rr = u_prim_permuted[j, ii, 1]
			v1_rr = u_prim_permuted[j, ii, 2]
			v2_rr = u_prim_permuted[j, ii, 3]
			theta_rr = u_prim_permuted[j, ii, 4]
			exner_rr = u_prim_permuted[j, ii, 5]
			phi_rr = u_prim_permuted[j, ii, 6]

			normal_direction_1 = 0.5 * (contravariant_vectors_x[j, i, 1] +
										contravariant_vectors_x[j, ii, 1])
			normal_direction_2 = 0.5 * (contravariant_vectors_x[j, i, 2] +
										contravariant_vectors_x[j, ii, 2])

			v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2
			v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2

			rho_avg = 0.5f0 * (rho_ll + rho_rr)
			v1_avg = 0.5f0 * (v1_ll + v1_rr)
			v2_avg = 0.5f0 * (v2_ll + v2_rr)
			theta_avg = 0.5f0 * (theta_ll + theta_rr)
			kin_avg = 0.5f0 *(v1_rr * v1_rr + v2_rr * v2_rr + v1_ll * v1_ll + v2_ll * v2_ll)
			# Calculate fluxes depending on normal_direction
			f1 = rho_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
			f2 = kin_avg * 0.5f0 * normal_direction_1
			f3 = kin_avg * 0.5f0 * normal_direction_2
			f4 = f1 * theta_avg
			gravity =  (phi_rr - phi_ll)
			jump_v1 = v1_rr - v1_ll
			jump_v2 = v2_rr - v2_ll
			theta_grad_exner = equations.c_p * theta_avg * (exner_rr - exner_ll)
			vorticity_x = v2_ll * jump_v1 * normal_direction_2 - v2_ll * jump_v2 * normal_direction_1			
			vorticity_y = v1_ll * jump_v2 * normal_direction_1 - v1_ll * jump_v1 * normal_direction_2
			g2 = 0.5f0 * vorticity_x + 0.5f0 * normal_direction_1 * (theta_grad_exner + gravity)
			g3 = 0.5f0 * vorticity_y + 0.5f0 * normal_direction_2 * (theta_grad_exner + gravity)

			# Add scaled fluxes to RHS
			factor_i = alpha * derivative_split[i, ii]
			du_permuted[j, i, 1] += factor_i * f1
			du_permuted[j, i, 2] += factor_i * (f2 + g2)
			du_permuted[j, i, 3] += factor_i * (f3 + g3)
			du_permuted[j, i, 4] += factor_i * f4

			factor_ii = alpha * derivative_split[ii, i]
			du_permuted[j, ii, 1] += factor_ii * f1
			du_permuted[j, ii, 2] += factor_ii * (f2 - g2)
			du_permuted[j, ii, 3] += factor_ii * (f3 - g3)
			du_permuted[j, ii, 4] += factor_ii * f4
		end
	end

	@turbo for v in eachvariable(equations),
		j in eachnode(dg),
		i in eachnode(dg)

		du[i, j, v] = du_permuted[j, i, v]
	end

	# y direction
	# We must also permute the contravariant vectors.
	contravariant_vectors_y = StrideArray{eltype(contravariant_vectors)}(undef,
		(StaticInt(nnodes(dg)),
			StaticInt(nnodes(dg)),
			StaticInt(ndims(mesh))))

	@turbo for k in eachnode(dg), j in eachnode(dg), i in eachnode(dg)
		contravariant_vectors_y[i, j, 1] = contravariant_vectors[1, 2, i, j, element]
		contravariant_vectors_y[i, j, 2] = contravariant_vectors[2, 2, i, j, element]
	end

	# The memory layout is already optimal for SIMD vectorization in this loop.
	for j in eachnode(dg), jj in (j+1):nnodes(dg)
		@turbo for i in eachnode(dg)
			rho_ll = u_prim[i, j, 1]
			v1_ll = u_prim[i, j, 2]
			v2_ll = u_prim[i, j, 3]
			theta_ll = u_prim[i, j, 4]
			exner_ll = u_prim[i, j, 5]
			phi_ll = u_prim[i, j, 6]

			rho_rr = u_prim[i, jj, 1]
			v1_rr = u_prim[i, jj, 2]
			v2_rr = u_prim[i, jj, 3]
			theta_rr = u_prim[i, jj, 4]
			exner_rr = u_prim[i, jj, 5]
			phi_rr = u_prim[i, jj, 6]

			normal_direction_1 = 0.5 * (contravariant_vectors_y[i, j, 1] +
										contravariant_vectors_y[i, jj, 1])
			normal_direction_2 = 0.5 * (contravariant_vectors_y[i, j, 2] +
										contravariant_vectors_y[i, jj, 2])


			v_dot_n_ll = v1_ll * normal_direction_1 + v2_ll * normal_direction_2
			v_dot_n_rr = v1_rr * normal_direction_1 + v2_rr * normal_direction_2

			rho_avg = 0.5f0 * (rho_ll + rho_rr)
			v1_avg = 0.5f0 * (v1_ll + v1_rr)
			v2_avg = 0.5f0 * (v2_ll + v2_rr)
			theta_avg = 0.5f0 * (theta_ll + theta_rr)
			kin_avg = 0.5f0 *(v1_rr * v1_rr + v2_rr * v2_rr + v1_ll * v1_ll + v2_ll * v2_ll)
			# Calculate fluxes depending on normal_direction
			f1 = rho_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
			f2 = kin_avg * 0.5f0 * normal_direction_1
			f3 = kin_avg * 0.5f0 * normal_direction_2
			f4 = f1 * theta_avg
			gravity =  (phi_rr - phi_ll)
			jump_v1 = v1_rr - v1_ll
			jump_v2 = v2_rr - v2_ll
			theta_grad_exner = equations.c_p * theta_avg * (exner_rr - exner_ll)
			vorticity_x = v2_ll * jump_v1 * normal_direction_2 - v2_ll * jump_v2 * normal_direction_1			
			vorticity_y = v1_ll * jump_v2 * normal_direction_1 - v1_ll * jump_v1 * normal_direction_2
			g2 = 0.5f0 * vorticity_x + 0.5f0 * normal_direction_1 * (theta_grad_exner + gravity)
			g3 = 0.5f0 * vorticity_y + 0.5f0 * normal_direction_2 * (theta_grad_exner + gravity)
			# Add scaled fluxes to RHS
			factor_j = alpha * derivative_split[j, jj]
			du[i, j, 1] += factor_j * f1
			du[i, j, 2] += factor_j * (f2 + g2)
			du[i, j, 3] += factor_j * (f3 + g3)
			du[i, j, 4] += factor_j * f4

			factor_jj = alpha * derivative_split[jj, j]
			du[i, jj, 1] += factor_jj * f1
			du[i, jj, 2] += factor_jj * (f2 - g2)
			du[i, jj, 3] += factor_jj * (f3 - g3)
			du[i, jj, 4] += factor_jj * f4
		end
	end

	# Finally, we add the temporary RHS computed here to the global RHS in the
	# given `element`.
	@turbo for v in eachvariable(equations),
		j in eachnode(dg),
		i in eachnode(dg)

		_du[v, i, j, element] += du[i, j, v]
	end
end
