# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

using Documenter, Mortar2D

makedocs(modules=[Mortar2D],
         format = Documenter.HTML(),
         sitename = "Mortar2D",
         pages = [
                  "Introduction" => "index.md",
                  "Theory" => "theory.md",
                  "API" => "api.md"
                 ])
