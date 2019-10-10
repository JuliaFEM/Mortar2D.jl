# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

using Test

test_files = String[]
push!(test_files, "test_calculate_normals.jl")
push!(test_files, "test_calculate_projections.jl")
push!(test_files, "test_calculate_segments.jl")
push!(test_files, "test_calculate_mortar_matrices.jl")
push!(test_files, "test_calculate_mortar_assembly.jl")

@testset "Mortar2D.jl" begin
    for fn in test_files
        include(fn)
    end
end
