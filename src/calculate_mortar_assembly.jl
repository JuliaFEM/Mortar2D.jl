# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

"""
    calculate_mortar_assembly(elements::Dict{Int, Vector{Int}},
                              element_types::Dict{Int, Symbol},
                              coords::Dict{Int, Vector{Float64}},
                              slave_element_ids::Vector{Int},
                              master_element_ids::Vector{Int})

Given data, calculate projection matrix `P`. This is the main function of
package.

Matrix ``P`` is defined as a projection between slave and master surface,
i.e.
```math
    D u_s = M u_m \\Rightarrow u_s = D^{-1} M u_m = P u_m.
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

    D = sparse(D_I, D_J, D_V)
    Df = cholfact(1/2*(D+D')[S,S])

    M_ = sparse(M_I, M_J, M_V)
    P =  Df \ M_[S,M]
    SparseArrays.droptol!(P, 1.0e-12)
    return P
end
