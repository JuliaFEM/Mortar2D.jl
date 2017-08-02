# Mortar2D.jl

[![Build Status](https://travis-ci.org/JuliaFEM/Mortar2D.jl.svg?branch=master)](https://travis-ci.org/JuliaFEM/Mortar2D.jl)[![Coverage Status](https://coveralls.io/repos/github/JuliaFEM/Mortar2D.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaFEM/Mortar2D.jl?branch=master)[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliafem.github.io/Mortar2D.jl/stable)[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliafem.github.io/Mortar2D.jl/latest)[![Issues](https://img.shields.io/github/issues/JuliaFEM/Mortar2D.jl.svg)](https://github.com/JuliaFEM/Mortar2D.jl/issues)

Mortar2D.jl is a Julia package to calculate discrete projections between
non-conforming finite element mesheds. The resulting "mortar matrices" can
be used to tie non-conforming finite element meshes together which are meshed
separately to construct bigger models.

Using mortar methods in mesh tie problems results variationally consistent
solution. Mathematically, goal is to solve mixed problem with primary field
variable and Lagrange multipliers, which have a physical meaning (e.g. contact
pressure if unknown field is displacement). The problem arising is a typical
saddle point problem with zeros on diagonal.

## Installing and testing package

Installing package goes same way like other packages in julia, i.e.
```julia
julia> Pkg.add("Mortar2D")
```

Testing package can be done using `Pkg.test`, i.e.
```julia
julia> Pkg.test("Mortar2D")
```

Probably the easiest way to test the functionality of package is to use [JuliaBox](https://juliabox.com/).
