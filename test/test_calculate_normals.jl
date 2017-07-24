# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

using Base.Test

using Mortar2D: calculate_normals

@testset "calculate normals for segments" begin
    elements = Dict(
        1 => [1, 2],
        2 => [2, 3])
    element_types = Dict(
        1 => :Seg2,
        2 => :Seg2)
    X = Dict(
        1 => [0.0, 0.0],
        2 => [1.0, 0.0],
        3 => [2.0, 0.0])
    n = calculate_normals(elements, element_types, X)
    @test n[1] == [0.0, 1.0]
    @test n[2] == [0.0, 1.0]
    @test n[3] == [0.0, 1.0]
end
