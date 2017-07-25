# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

"""
Given slave surface elements, master surface elements, nodal coordinates and
normal direction on nodes of slave surface elements, calculate contact
segments.

Returns a dictionary, where key is slave element id and values is a list of
master elements giving contribution to that slave elements and xi-coordinates
of slave side element.
"""
function calculate_segments(slave_elements, master_elements, element_types, coords, normals)
    S = Dict()
    for (sid, scon) in slave_elements
        @assert element_types[sid] == :Seg2
        xs1 = coords[scon[1]]
        xs2 = coords[scon[2]]
        ns1 = normals[scon[1]]
        ns2 = normals[scon[2]]
        for (mid, mcon) in master_elements
            @assert element_types[mid] == :Seg2
            xm1 = coords[mcon[1]]
            xm2 = coords[mcon[2]]
            xi1a = project_from_master_to_slave(Val{:Seg2}, xm1, xs1, xs2, ns1, ns2)
            xi1b = project_from_master_to_slave(Val{:Seg2}, xm2, xs1, xs2, ns1, ns2)
            xi1 = clamp.([xi1a; xi1b], -1.0, 1.0)
            l = 1/2*abs(xi1[2]-xi1[1])
            if isapprox(l, 0.0)
                # no contribution from master so slave, we are done
                continue
            end
            if !haskey(S, sid)
                S[sid] = []
            end
            push!(S[sid], (mid, xi1))
        end
    end
    return S
end
