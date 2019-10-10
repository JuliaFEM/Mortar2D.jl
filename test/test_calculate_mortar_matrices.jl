# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

using Test
using Mortar2D: calculate_mortar_matrices

@testset "calculate mortar matrices" begin
    elements = Dict(
        1 => [1, 2],
        2 => [3, 4])
    element_types = Dict(
        1 => :Seg2,
        2 => :Seg2)
    coords = Dict(
        1 => [1.0, 2.0],
        2 => [3.0, 2.0],
        3 => [2.0, 2.0],
        4 => [0.0, 2.0])
    normals = Dict(
        1 => [0.0, -1.0],
        2 => [0.0, -1.0])
    segmentation = Dict(1 => [(2, [-1.0, 0.0])])
    De, Me = calculate_mortar_matrices(1, elements, element_types,
                                       coords, normals, segmentation)
    De_expected = 1/24*[14 4; 4 2]
    Me_expected = 1/24*[13 5; 5 1]
    @test isapprox(De, De_expected)
    @test isapprox(Me[2], Me_expected)
end
