# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

using Test

@testset "Mortar2D.jl" begin
    @testset "test_calculate_normals.jl" begin include("test_calculate_normals.jl") end
    @testset "test_calculate_projections.jl" begin include("test_calculate_projections.jl") end
    @testset "test_calculate_segments.jl" begin include("test_calculate_segments.jl") end
    @testset "test_calculate_mortar_matrices.jl" begin include("test_calculate_mortar_matrices.jl") end
    @testset "test_calculate_mortar_assembly.jl" begin include("test_calculate_mortar_assembly.jl") end
end
