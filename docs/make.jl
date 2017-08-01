# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

using Documenter, Mortar2D

installed_packages = Pkg.installed()
if !haskey(installed_packages, "PyPlot")
    Pkg.add("PyPlot")
end

makedocs(modules=[Mortar2D],
         format = :html,
         sitename = "Mortar2D",
         pages = [
                  "Introduction" => "index.md",
                  "Theory" => "theory.md",
                  "API" => "api.md"
                 ])
