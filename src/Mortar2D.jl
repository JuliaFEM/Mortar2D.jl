# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

module Mortar2D
include("calculate_normals.jl")
include("project_nodes.jl")
include("calculate_segments.jl")
end
