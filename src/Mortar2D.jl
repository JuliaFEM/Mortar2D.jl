# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

module Mortar2D
include("calculate_normals.jl")
include("calculate_projections.jl")
include("calculate_segments.jl")
include("calculate_mortar_matrices.jl")
include("calculate_assembly.jl")
end
