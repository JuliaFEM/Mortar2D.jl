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

Calculate mortar matrices for simple problem in README.md

```jldoctest ex1
Xs = Dict(1 => [0.0, 1.0], 2 => [5/4, 1.0], 3 => [2.0, 1.0])
Xm = Dict(4 => [0.0, 1.0], 5 => [1.0, 1.0], 6 => [2.0, 1.0])
coords = merge(Xm , Xs)
Es = Dict(1 => [1, 2], 2 => [2, 3])
Em = Dict(3 => [4, 5], 4 => [5, 6])
elements = merge(Es, Em)
element_types = Dict(1 => :Seg2, 2 => :Seg2, 3 => :Seg2, 4 => :Seg2)
slave_element_ids = [1, 2]
master_element_ids = [3, 4]
s, m, D, M = calculate_mortar_assembly(elements, element_types, coords,
                                       slave_element_ids, master_element_ids)

# output

([1, 2, 3], [4, 5, 6],
  [1, 1]  =  0.416667
  [2, 1]  =  0.208333
  [1, 2]  =  0.208333
  [2, 2]  =  0.666667
  [3, 2]  =  0.125
  [2, 3]  =  0.125
  [3, 3]  =  0.25,
  [1, 4]  =  0.366667
  [2, 4]  =  0.133333
  [1, 5]  =  0.25625
  [2, 5]  =  0.65
  [3, 5]  =  0.09375
  [1, 6]  =  0.00208333
  [2, 6]  =  0.216667
  [3, 6]  =  0.28125)

```

`s` and `m` contains slave and master dofs:
```jldoctest ex1
julia> s, m
([1, 2, 3], [4, 5, 6])
```

`D` is slave side mortar matrix:
```jldoctest ex1
julia> full(D[s,s])
3×3 Array{Float64,2}:
 0.416667  0.208333  0.0
 0.208333  0.666667  0.125
 0.0       0.125     0.25
```

`M` is master side mortar matrix:
```jldoctest ex1
julia> full(M[s,m])
3×3 Array{Float64,2}:
 0.366667  0.25625  0.00208333
 0.133333  0.65     0.216667
 0.0       0.09375  0.28125
```

"""
function calculate_mortar_assembly(elements, element_types, coords,
                                   slave_element_ids, master_element_ids)
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
