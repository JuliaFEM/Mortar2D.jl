# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

function assemble!(problem::Problem{Mortar}, time::Float64, ::Type{Val{1}}, ::Type{Val{false}})

    props = problem.properties
    field_dim = get_unknown_field_dimension(problem)
    field_name = get_parent_field_name(problem)
    slave_elements = get_slave_elements(problem)

    # 1. calculate nodal normals and tangents for slave element nodes j ∈ S
    normals, tangents = calculate_normals(slave_elements, time, Val{1};
                                          rotate_normals=props.rotate_normals)
    update!(slave_elements, "normal", time => normals)
    update!(slave_elements, "tangent", time => tangents)

    # 2. loop all slave elements
    for slave_element in slave_elements

        nsl = length(slave_element)
        X1 = slave_element("geometry", time)
        n1 = slave_element("normal", time)

        # 3. loop all master elements
        for master_element in slave_element("master elements", time)

            nm = length(master_element)
            X2 = master_element("geometry", time)

            # 3.1 calculate segmentation
            xi1a = project_from_master_to_slave(slave_element, X2[1], time)
            xi1b = project_from_master_to_slave(slave_element, X2[2], time)
            xi1 = clamp.([xi1a; xi1b], -1.0, 1.0)
            l = 1/2*abs(xi1[2]-xi1[1])
            isapprox(l, 0.0) && continue # no contribution in this master element

            # 3.2. bi-orthogonal basis
            De = zeros(nsl, nsl)
            Me = zeros(nsl, nsl)
            Ae = zeros(nsl, nsl)
            if props.dual_basis
                for ip in get_integration_points(slave_element, 3)
                    detJ = slave_element(ip, time, Val{:detJ})
                    w = ip.weight*detJ*l
                    xi = ip.coords[1]
                    xi_s = dot([1/2*(1-xi); 1/2*(1+xi)], xi1)
                    N1 = vec(get_basis(slave_element, xi_s, time))
                    De += w*diagm(N1)
                    Me += w*N1*N1'
                end
                Ae = De*inv(Me)
            else
                Ae = eye(nsl)
            end

            # 3.3. loop integration points of one integration segment and calculate
            # local mortar matrices
            fill!(De, 0.0)
            fill!(Me, 0.0)
            ge = zeros(field_dim*nsl)
            for ip in get_integration_points(slave_element, 2)
                detJ = slave_element(ip, time, Val{:detJ})
                w = ip.weight*detJ*l
                xi = ip.coords[1]
                xi_s = dot([1/2*(1-xi); 1/2*(1+xi)], xi1)
                N1 = vec(get_basis(slave_element, xi_s, time))
                Phi = Ae*N1
                # project gauss point from slave element to master element in direction n_s
                X_s = N1*X1 # coordinate in gauss point
                n_s = N1*n1 # normal direction in gauss point
                xi_m = project_from_slave_to_master(master_element, X_s, n_s, time)
                N2 = vec(get_basis(master_element, xi_m, time))
                X_m = N2*X2
                De += w*Phi*N1'
                Me += w*Phi*N2'
                if props.adjust
                    haskey(slave_element, "displacement") || continue
                    haskey(master_element, "displacement") || continue
                    norm(mean(X1) - X2[1]) / norm(X1[2] - X1[1]) < props.distval || continue
                    norm(mean(X1) - X2[2]) / norm(X1[2] - X1[1]) < props.distval || continue
                    u1 = slave_element("displacement", time)
                    u2 = master_element("displacement", time)
                    x_s = X_s + N1*u1
                    x_m = X_m + N2*u2
                    ge += w*vec((x_m-x_s)*Phi')
                end
            end

            # add contribution to contact virtual work
            sdofs = get_gdofs(problem, slave_element)
            mdofs = get_gdofs(problem, master_element)

            for i=1:field_dim
                lsdofs = sdofs[i:field_dim:end]
                lmdofs = mdofs[i:field_dim:end]
                add!(problem.assembly.C1, lsdofs, lsdofs, De)
                add!(problem.assembly.C1, lsdofs, lmdofs, -Me)
                add!(problem.assembly.C2, lsdofs, lsdofs, De)
                add!(problem.assembly.C2, lsdofs, lmdofs, -Me)
            end
            add!(problem.assembly.g, sdofs, ge)

        end # master elements done

    end # slave elements done, contact virtual work ready

end
