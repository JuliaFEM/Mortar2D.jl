# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

const MortarSegmentation = Dict{Int, Vector{Tuple{Int, Vector{Float64}}}}

"""
    calculate_segments(slave_element_ids::Vector{Int},
                       master_element_ids::Vector{Int},
                       elements::Dict{Int, Vector{Int}},
                       element_types::Dict{Int, Symbol},
                       coords::Dict{Int, Vector{Float64}},
                       normals::Dict{Int, Vector{Float64}})

Given slave surface elements, master surface elements, nodal coordinates and
normal direction on nodes of slave surface elements, calculate contact
segments.

Return type is a dictionary, where key is slave element id and values is a
list of master elements giving contribution to that slave elements and
xi-coordinates of slave side element.

# Example
```jldoctest
elements = Dict(1 => [1, 2], 2 => [3, 4])
element_types = Dict(1 => :Seg2, 2 => :Seg2)
coords = Dict(
    1 => [1.0, 2.0],
    2 => [3.0, 2.0],
    3 => [0.0, 2.0],
    4 => [2.0, 2.0])
normals = Dict(
    1 => [0.0, -1.0],
    2 => [0.0, -1.0])
slave_ids = [1]
master_ids = [2]
segments = calculate_segments(slave_ids, master_ids, elements,
                              element_types, coords, normals)

# output

Dict{Int64,Array{Tuple{Int64,Array{Float64,1}},1}} with 1 entry:
  1 => Tuple{Int64,Array{Float64,1}}[(2, [-1.0, -0.0])]

```

Here, output result means that slave element #1 has segment with master
element(s) #2 with dimensionless slave element coordinate xi = [-1, 0].
That is, the start and end point of projection in physical coordinate
system is:

    x_start = 1/2*(1-xi[1])*xs1 + 1/2*(1+xi[1])*xs2
    x_stop = 1/2*(1-xi[2])*xs1 + 1/2*(1+xi[2])*xs2

"""
function calculate_segments(slave_element_ids::Vector{Int},
                            master_element_ids::Vector{Int},
                            elements::Dict{Int, Vector{Int}},
                            element_types::Dict{Int, Symbol},
                            coords::Dict{Int, Vector{Float64}},
                            normals::Dict{Int, Vector{Float64}})
    S = MortarSegmentation()
    for sid in slave_element_ids
        @assert element_types[sid] == :Seg2
        haskey(S, sid) || (S[sid] = [])
        scon = elements[sid]
        xs1 = coords[scon[1]]
        xs2 = coords[scon[2]]
        ns1 = normals[scon[1]]
        ns2 = normals[scon[2]]
        for mid in master_element_ids
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
