# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

using Test
using Mortar2D: calculate_segments

# calculate segmentation
elements = Dict(
    1 => [1, 2],
    2 => [3, 4])
element_types = Dict(
    1 => :Seg2,
    2 => :Seg2)
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
segments = calculate_segments(
    slave_ids,
    master_ids,
    elements,
    element_types,
    coords,
    normals)
@test length(segments) == 1
@test length(segments[1]) == 1
elid, xi = first(segments[1])
@test elid == 2
@test isapprox(xi, [-1.0, 0.0])
