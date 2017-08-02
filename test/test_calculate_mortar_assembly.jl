# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

using Base.Test
using Mortar2D: calculate_mortar_assembly

@testset "calculate mortar assembly" begin
    ns = 10 # number of slaves in interface
    nm = 9 # number of masters in interface

    Xs = Dict(i => [x, 1/2] for (i, x) in enumerate(linspace(0, 1, ns)))
    Es = Dict(i => [i, i+1] for i=1:(ns-1))
    Xm = Dict(ns+i => [x, 1/2] for (i, x) in enumerate(linspace(0, 1, nm)))
    Em = Dict(i => [i, i+1]+1 for i=ns:ns+nm-2)
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
    um = linspace(0, 1, nm)
    us = D[s,s] \ (M[s,m]*um)
    us_expected = linspace(0, 1, ns)
    @test isapprox(us, us_expected)
end
