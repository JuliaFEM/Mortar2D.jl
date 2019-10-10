# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

"""
    calculate_normals(elements::Dict{Int, Vector{Int}},
                      element_types::Dict{Int, Symbol},
                      X::Dict{Int, Vector{Float64})

Given elements, element types and node locations, calculate nodal normals by
first calculating normal directions for each element and then averaging them
in nodes. As a result we get unique normal direction defined to each node.

# Notes
Only linear elements supported.

# Example

```jldoctest
X = Dict(1 => [7.0, 7.0], 2 => [4.0, 3.0], 3 => [0.0, 0.0])
elements = Dict(1 => [1, 2], 2 => [2, 3])
element_types = Dict(1 => :Seg2, 2 => :Seg2)
normals = calculate_normals(elements, element_types, X)

# output

Dict{Int64,Array{Float64,1}} with 3 entries:
  2 => [0.707107, -0.707107]
  3 => [0.6, -0.8]
  1 => [0.8, -0.6]

```

"""
function calculate_normals(elements, element_types, X)
    normals = empty(X)
    for (elid, elcon) in elements
        @assert element_types[elid] == :Seg2
        d = X[elcon[2]] - X[elcon[1]]
        n = [-d[2], d[1]]
        n /= norm(n)
        for nid in elcon
            haskey(normals, nid) || (normals[nid] = zeros(2))
            normals[nid] += n
        end
    end
    for (nid, n) in normals
        normals[nid] /= norm(normals[nid])
    end
    return normals
end
