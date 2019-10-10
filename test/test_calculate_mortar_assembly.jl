# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

using Test
using Mortar2D: calculate_mortar_assembly

# calculate mortar assembly
ns = 10 # number of slaves in interface
nm = 9 # number of masters in interface

Xs = Dict(i => [x, 1/2] for (i, x) in enumerate(range(0.0, stop=1.0, length=ns)))
Es = Dict(i => [i, i+1] for i=1:(ns-1))
Xm = Dict(ns+i => [x, 1/2] for (i, x) in enumerate(range(0.0, stop=1.0, length=nm)))
Em = Dict(i => [i+1, i+2] for i=ns:ns+nm-2)
element_types = Dict(i => :Seg2 for i=1:(ns+nm))

slave_element_ids = collect(keys(Es))
master_element_ids = collect(keys(Em))
elements = merge(Es, Em)
coords = merge(Xs, Xm)

s, m, D, M = calculate_mortar_assembly(elements,
                                       element_types,
                                       coords,
                                       slave_element_ids,
                                       master_element_ids)

# test case 1, constant pressure load
um = ones(nm)
us = D[s,s] \ (M[s,m]*um)
us_expected = ones(ns)
@test isapprox(us, us_expected)

# test case 2, linear load
um = range(0.0, stop=1.0, length=nm)
us = D[s,s] \ (M[s,m]*um)
us_expected = range(0.0, stop=1.0, length=ns)
@test isapprox(us, us_expected)
