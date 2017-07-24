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
        1 => [7.0, 7.0],
        2 => [4.0, 3.0],
        3 => [0.0, 0.0])
    n = calculate_normals(elements, element_types, X)
    @test isapprox(n[1], 1/5*[4, -3])
    @test isapprox(n[2], sqrt(2)/2*[1, -1])
    @test isapprox(n[3], 1/5*[3, -4])
end
