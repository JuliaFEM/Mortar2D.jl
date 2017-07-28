# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

"""
    calculate_normals(elements::Dict{Int, Vector{Int}}, element_types::Dict{Int, Symbol}, X::Dict{Int, Vector{Float64})

Given elements, element types and node locations, calculate nodal normals by
first calculating normal directions for each element and then averaging them
in nodes. As a result we get unique normal direction defined to each node.

# References

- Yang2005
"""
function calculate_normals{T<:Integer,P<:AbstractFloat}(elements::Dict{T, Vector{T}},
                                                        element_types::Dict{T, Symbol},
                                                        X::Dict{T, Vector{P}})
    normals = Dict{T, Vector{P}}()
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
