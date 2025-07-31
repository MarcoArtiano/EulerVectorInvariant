using Trixi: @threaded, False, VolumeIntegralWeakForm, DGSEM, flux, eachelement, eachnode, get_node_vars, multiply_add_to_node_vars!, get_contravariant_vector
import Trixi: calc_volume_integral!

@muladd begin

	function calc_volume_integral!(du, u,
		mesh::Union{TreeMesh{2}, StructuredMesh{2},
			StructuredMeshView{2}, UnstructuredMesh2D,
			P4estMesh{2}, P4estMeshView{2},
			T8codeMesh{2}},
		nonconservative_terms, equations::CompressibleEulerVectorInvariantEquations2D,
		volume_integral::VolumeIntegralWeakForm,
		dg::DGSEM, cache)

		@threaded for element in eachelement(dg, cache)
			strong_form_kernel!(du, u, element, mesh,
				nonconservative_terms, equations,
				dg, cache)
		end

		return nothing
	end

	@inline function get_covariant_vector(index, contravariant_vectors, i, j, element)
		sign = 3 - 2 * index
		return SVector(-sign * contravariant_vectors[index, 2, i, j, element],
			sign * contravariant_vectors[index, 1, i, j, element])
	end


	@inline function strong_form_kernel!(du, u,
		element,
		mesh::Union{StructuredMesh{2}, StructuredMeshView{2},
			UnstructuredMesh2D, P4estMesh{2},
			P4estMeshView{2}, T8codeMesh{2}},
		nonconservative_terms::False, equations::CompressibleEulerVectorInvariantEquations2D,
		dg::DGSEM, cache, alpha = true)
		# true * [some floating point value] == [exactly the same floating point value]
		# This can (hopefully) be optimized away due to constant propagation.
		@unpack derivative_split, weights = dg.basis
		@unpack contravariant_vectors = cache.elements
		## Correction of the derivative operator for the strong form
		derivative_split[1, 1] -= 1 / weights[1]
		derivative_split[end, end] += 1 / weights[end]
		@. derivative_split = derivative_split * 0.5
		derivative_split[1, 1] += 1 / weights[1]
		derivative_split[end, end] -= 1 / weights[end]

		for j in eachnode(dg), i in eachnode(dg)
			u_node = get_node_vars(u, equations, dg, i, j, element)
			## rho, v1, v2, rho_theta
			exner = exner_pressure(u_node, equations)
			flux1 = flux(u_node, 1, equations) #rho_v1, 0 , 0, rho_theta_v1
			flux2 = flux(u_node, 2, equations) #rho_v2, 0 , 0, rho_theta_v2
			rho = u_node[1]
			theta = u_node[4] / u_node[1]

			# u1, u2
			## X Direction
			# Compute the contravariant flux by taking the scalar product of the
			# first contravariant vector Ja^1 and the flux vector
			Ja11, Ja12 = get_contravariant_vector(1, contravariant_vectors, i, j, element)
			contravariant_flux1 = Ja11 * flux1 + Ja12 * flux2

			## Y Direction
			# Compute the contravariant flux by taking the scalar product of the
			# second contravariant vector Ja^2 and the flux vector
			Ja21, Ja22 = get_contravariant_vector(2, contravariant_vectors, i, j, element)
			contravariant_flux2 = Ja21 * flux1 + Ja22 * flux2

			a11, a12 = get_covariant_vector(2, contravariant_vectors, i, j, element)
			a21, a22 = get_covariant_vector(1, contravariant_vectors, i, j, element)

			u_covariant_1 = a11 * u_node[2] + a12 * u_node[3]
			u_covariant_2 = a21 * u_node[2] + a22 * u_node[3]
            		u_covariant_1 = u_node[2]
            		u_covariant_2 = u_node[3]
			normu = contravariant_flux1[1] / rho * u_covariant_1 + contravariant_flux2[1] / rho * u_covariant_2
			
			for ii in eachnode(dg)
				## Density
				du[1, ii, j, element] = du[1, ii, j, element] + derivative_split[ii, i] * contravariant_flux1[1]

				## Potential Temperature: conservative
				du[4, ii, j, element] = du[4, ii, j, element] + 0.5 * derivative_split[ii, i] * contravariant_flux1[4]

				## Potential Temperature: non - conservative
	#			du[4, ii, j, element] = du[4, ii, j, element] + 0.5 * theta * derivative_split[ii, i] * contravariant_flux1[1] + 0.5 * contravariant_flux1[1] * derivative_split[ii, i] * theta

				## Momentum: conservative
				du[2, ii, j, element] = du[2, ii, j, element] + 0.5 * derivative_split[ii, i] * normu
				
			#	du[2, ii, j, element] = du[2, ii, j, element] - contravariant_flux2[1] / rho * derivative_split[ii, i] * u_covariant_2 + 0.5 * derivative_split[ii, i] * normu + (theta * derivative_split[ii, i] * exner) * Ja11

			#	du[3, ii, j, element] = du[3, ii, j, element] + contravariant_flux1[1] / rho * derivative_split[ii, i] * u_covariant_1
			end

			for jj in eachnode(dg)
				## Density
				du[1, i, jj, element] = du[1, i, jj, element] + derivative_split[jj, j] * contravariant_flux2[1]

				## Potential Temperature: conservative
				du[4, i, jj, element] = du[4, i, jj, element] + 0.5 * derivative_split[jj, j] * contravariant_flux2[4]

				## Potential Temperature: non - conservative
	#			du[4, i, jj, element] = du[4, i, jj, element] + 0.5 * theta * derivative_split[jj, j] * contravariant_flux2[1] + 0.5 * contravariant_flux2[1] * derivative_split[jj, j] * theta

				## Momentum: conservative
				du[3, i, jj, element] = du[3, i, jj, element] + 0.5 * derivative_split[jj, j] * normu
			#	du[2, i, jj, element] = du[2, i, jj, element] + contravariant_flux2[1] / rho * derivative_split[jj, j] * u_covariant_2

			#	du[3, i, jj, element] = du[3, i, jj, element] - contravariant_flux1[1] / rho * derivative_split[jj, j] * u_covariant_1 + 0.5 * derivative_split[jj, j] * normu + (theta * derivative_split[jj, j] * exner) * Ja22
			end
		end
	## Compute the non-conservative terms...
		for j in eachnode(dg), i in eachnode(dg)
		## X Direction
		tmp1 = 0.0
		tmp2 = 0.0
		tmp3 = 0.0
		tmp4 = 0.0
		tmp5 = 0.0
		for k in eachnode(dg)
		u_node = get_node_vars(u, equations, dg, k, j, element)
		flux1 = flux(u_node, 1, equations) #rho_v1, 0 , 0, rho_theta_v1
		flux2 = flux(u_node, 2, equations) #rho_v2, 0 , 0, rho_theta_v2
		exner = exner_pressure(u_node, equations)
		theta = u_node[4] / u_node[1]
		Ja11, Ja12 = get_contravariant_vector(1, contravariant_vectors, k, j, element)
		contravariant_flux1 = Ja11 * flux1 + Ja12 * flux2
		Ja21, Ja22 = get_contravariant_vector(2, contravariant_vectors, k, j, element)
		contravariant_flux2 = Ja21 * flux1 + Ja22 * flux2
		tmp1 +=  derivative_split[i, k] * u_node[1]
		tmp2 +=  derivative_split[i, k] * u_node[2]
		tmp3 +=  derivative_split[i, k] * exner
		tmp4 +=  derivative_split[i, k] * theta
		tmp5 +=  derivative_split[i, k] * contravariant_flux1[1]
		end

		u_node = get_node_vars(u, equations, dg, i, j, element)
		flux1 = flux(u_node, 1, equations) #rho_v1, 0 , 0, rho_theta_v1
		flux2 = flux(u_node, 2, equations) #rho_v2, 0 , 0, rho_theta_v2
		rho = u_node[1]
		theta = u_node[4] / u_node[1]
		Ja11, Ja12 = get_contravariant_vector(1, contravariant_vectors, i, j, element)
		contravariant_flux1 = Ja11 * flux1 + Ja12 * flux2
		Ja21, Ja22 = get_contravariant_vector(2, contravariant_vectors, i, j, element)
		contravariant_flux2 = Ja21 * flux1 + Ja22 * flux2
		du[2,i,j,element] += -contravariant_flux2[1]/rho * tmp2 + theta * tmp3 * Ja11
		du[3,i,j,element] += contravariant_flux1[1]/rho * tmp2
		du[4,i,j,element] += contravariant_flux1[1] * tmp4 * 0.5 + theta * tmp5 * 0.5	
		## Y Direction
		tmp1 = 0.0
		tmp2 = 0.0
		tmp3 = 0.0
		tmp4 = 0.0
		tmp5 = 0.0
		for k in eachnode(dg)
		u_node = get_node_vars(u, equations, dg, i, k, element)
		exner = exner_pressure(u_node, equations)
		theta = u_node[4] / u_node[1]
		Ja11, Ja12 = get_contravariant_vector(1, contravariant_vectors, i, k, element)
		contravariant_flux1 = Ja11 * flux1 + Ja12 * flux2
		Ja21, Ja22 = get_contravariant_vector(2, contravariant_vectors, i, k, element)
		contravariant_flux2 = Ja21 * flux1 + Ja22 * flux2
		tmp1 += derivative_split[j, k] * u_node[1]
		tmp2 += derivative_split[j, k] * u_node[2]
		tmp3 += derivative_split[j, k] * exner
		tmp4 += derivative_split[j, k] * theta
		tmp5 += derivative_split[j, k] * contravariant_flux2[1]
		end

		u_node = get_node_vars(u, equations, dg, i, j, element)
		flux1 = flux(u_node, 1, equations) #rho_v1, 0 , 0, rho_theta_v1
		flux2 = flux(u_node, 2, equations) #rho_v2, 0 , 0, rho_theta_v2
		rho = u_node[1]
		theta = u_node[4] / u_node[1]
		Ja11, Ja12 = get_contravariant_vector(1, contravariant_vectors, i, j, element)
		contravariant_flux1 = Ja11 * flux1 + Ja12 * flux2
		Ja21, Ja22 = get_contravariant_vector(2, contravariant_vectors, i, j, element)
		contravariant_flux2 = Ja21 * flux1 + Ja22 * flux2
		du[2,i,j,element] += contravariant_flux2[1]/rho * tmp1
		du[3,i,j,element] += -contravariant_flux1[1]/rho * tmp1 + theta * tmp3 * Ja22
		du[4,i,j,element] += contravariant_flux2[1] * tmp4 * 0.5 + theta * tmp5 * 0.5
		end
		
		return nothing
	end

end
