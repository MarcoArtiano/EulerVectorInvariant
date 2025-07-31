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
    sign = 3 - 2*index
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
    @unpack derivative_dhat = dg.basis
    @unpack contravariant_vectors = cache.elements
    ## Correction of the derivative operator for the strong form
    derivative_dhat[1, 1] -= 1 / weights[1]
	derivative_dhat[end, end] += 1 / weights[end]
    @. derivative_dhat = derivative_dhat*0.5
    derivative_dhat[1, 1] += 1 / weights[1]
	derivative_dhat[end, end] -= 1 / weights[end]

    for j in eachnode(dg), i in eachnode(dg)
        u_node = get_node_vars(u, equations, dg, i, j, element)
        ## rho, v1, v2, rho_theta
        exner = exner_pressure(u_node, equations)
        flux1 = flux(u_node, 1, equations) #rho_v1, 0 , 0, rho_theta_v1
        flux2 = flux(u_node, 2, equations) #rho_v2, 0 , 0, rho_theta_v2
        rho = u_node[1]
        theta = u_node[4]/u_node[1]
        
        # u1, u2
        ## X Direction
        # Compute the contravariant flux by taking the scalar product of the
        # first contravariant vector Ja^1 and the flux vector
        Ja11, Ja12 = get_contravariant_vector(1, contravariant_vectors, i, j, element)
        contravariant_flux1 = Ja11 * flux1 + Ja12 * flux2 
        
        a11, a12 = get_covariant_vector(1, contravariant_vectors, i, j, element)
        a21, a22 = get_covariant_vector(2, contravariant_vectors, i, j, element)

        u_covariant_1 = a11 * u_node[2] + a12 * u_node[3]
        u_covariant_2 = a21 * u_node[2] + a22 * u_node[3]

        for ii in eachnode(dg)
            ## Density
            du[1, ii, j, element] = du[1, ii, j, element] + derivative_dhat[ii, i] * contravariant_flux1[1]

            ## Potential Temperature: conservative
            du[4, ii, j, element] = du[4, ii, j, element] + 0.5 * derivative_dhat[ii, i] * contravariant_flux1[4]

            ## Potential Temperature: non - conservative
            du[4, ii, j, element] = du[4, ii, j, element] + 0.5 * theta * derivative_dhat[ii, i] * contravariant_flux1[1] + 0.5 * contravariant_flux1[1] * derivative_dhat[ii, i] * theta

            normu = contravariant_flux1[1]/rho * u_covariant_1 + contravariant_flux1[2]/rho * u_covariant_2
            ## Momentum: conservative
            du[2, ii, j, element] = du[2, ii, j, element] - contravariant_flux1[2]/rho * derivative_dhat[ii, i] * u_covariant_2 + 0.5 * derivative_dhat[ii, i] * normu + theta * derivative_dhat[ii, i] * exner

            du[3, ii, j, element] = du[3, ii, j, element] + contravariant_flux1[1]/rho * derivative_dhat[ii, i] * u_covariant_2

            multiply_add_to_node_vars!(du, alpha * derivative_dhat[ii, i],
                                       contravariant_flux1, equations, dg,
                                       ii, j, element)
        end

        ## Y Direction

        # Compute the contravariant flux by taking the scalar product of the
        # second contravariant vector Ja^2 and the flux vector
        Ja21, Ja22 = get_contravariant_vector(2, contravariant_vectors, i, j, element)
        contravariant_flux2 = Ja21 * flux1 + Ja22 * flux2
        normu = contravariant_flux1[1]/rho * u_covariant_1 + contravariant_flux2[1]/rho * u_covariant_2

        for jj in eachnode(dg)
            ## Density
            du[1, i, jj, element] = du[1, i, jj, element] + derivative_dhat[jj, j] * contravariant_flux2[1]

            ## Potential Temperature: conservative
            du[4, i, jj, element] = du[4, i, jj, element] + 0.5 * derivative_dhat[jj, j] * contravariant_flux2[4]

            ## Potential Temperature: non - conservative
            du[4, i, jj, element] = du[4, i, jj, element] + 0.5 * theta * derivative_dhat[jj, j] * contravariant_flux2[1] + 0.5 * contravariant_flux2[1] * derivative_dhat[jj, j] * theta

            ## Momentum:
            du[2, i, jj, element] = du[2, i, jj, element] + contravariant_flux2[1]/rho * derivative_dhat[jj, j] * u_covariant_1 

            du[3, i, jj, element] = du[3, i, jj, element] - contravariant_flux1[1]/rho * derivative_dhat[jj, j] * u_covariant_1 + 0.5 * derivative_dhat[jj, j] * normu + theta * derivative_dhat[jj, j] * exner
        end
    end

    return nothing
end