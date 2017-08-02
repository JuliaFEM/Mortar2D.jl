# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

"""
    calculate_mortar_assembly(elements::Dict{Int, Vector{Int}},
                              element_types::Dict{Int, Symbol},
                              coords::Dict{Int, Vector{Float64}},
                              slave_element_ids::Vector{Int},
                              master_element_ids::Vector{Int})

Given data, calculate projection matrices `D` and `M`. This is the main
function of package. Relation between matrices is ``D u_s = M u_m``, where
``u_s`` is slave nodes and ``u_m`` master nodes.

# Example

```julia
s, m, D, M = calculate_mortar_assembly(elements, element_types, coords,
                                       slave_element_ids, master_element_ids)
```

"""
function calculate_mortar_assembly(elements::Dict{Int, Vector{Int}},
                                   element_types::Dict{Int, Symbol},
                                   coords::Dict{Int, Vector{Float64}},
                                   slave_element_ids::Vector{Int},
                                   master_element_ids::Vector{Int})
    S = Int[]
    M = Int[]
    for sid in slave_element_ids
        push!(S, elements[sid]...)
    end
    for mid in master_element_ids
        push!(M, elements[mid]...)
    end
    S = sort(unique(S))
    M = sort(unique(M))
    ns = length(S)
    nm = length(M)
    normals = calculate_normals(elements, element_types, coords)
    segmentation = calculate_segments(slave_element_ids, master_element_ids,
                                      elements, element_types, coords, normals)
    D_I = Int[]
    D_J = Int[]
    D_V = Float64[]
    M_I = Int[]
    M_J = Int[]
    M_V = Float64[]
    for sid in slave_element_ids
        mids = [mid for (mid, xi) in segmentation[sid]]
        De, Me = calculate_mortar_matrices(sid, elements, element_types,
                                           coords, normals, segmentation)
        for (i,I) in enumerate(elements[sid])
            for (j,J) in enumerate(elements[sid])
                push!(D_I, I)
                push!(D_J, J)
                push!(D_V, De[i,j])
            end
        end
        for mid in mids
            for (i,I) in enumerate(elements[sid])
                for (j,J) in enumerate(elements[mid])
                   push!(M_I, I)
                   push!(M_J, J)
                   push!(M_V, Me[mid][i,j])
                end
            end
        end
    end

    # return global matrices + slave and master dofs
    Dg = sparse(D_I, D_J, D_V)
    Mg = sparse(M_I, M_J, M_V)
    return S, M, Dg, Mg
end
