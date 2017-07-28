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
function calculate_segments(slave_ids, master_ids, elements, element_types, coords, normals)
    S = Dict()
    for sid in slave_ids
        @assert element_types[sid] == :Seg2
        haskey(S, sid) || (S[sid] = [])
        scon = elements[sid]
        xs1 = coords[scon[1]]
        xs2 = coords[scon[2]]
        ns1 = normals[scon[1]]
        ns2 = normals[scon[2]]
        for mid in master_ids
            @assert element_types[mid] == :Seg2
            mcon = elements[mid]
            xm1 = coords[mcon[1]]
            xm2 = coords[mcon[2]]
            # first project from slave to master, to find out
            # are we completely outside of domain
            xi2a = project_from_slave_to_master(Val{:Seg2}, xs1, ns1, xm1, xm2)
            xi2b = project_from_slave_to_master(Val{:Seg2}, xs2, ns2, xm1, xm2)
            xi2a >  1.0 && xi2b >  1.0 && continue
            xi2a < -1.0 && xi2b < -1.0 && continue
            xi1a = project_from_master_to_slave(Val{:Seg2}, xm1, xs1, xs2, ns1, ns2)
            xi1b = project_from_master_to_slave(Val{:Seg2}, xm2, xs1, xs2, ns1, ns2)
            xi1 = clamp.([xi1a; xi1b], -1.0, 1.0)
            isapprox(abs(xi1[2]-xi1[1]), 0.0) && continue
            push!(S[sid], (mid, xi1))
        end
    end
    return S
end
