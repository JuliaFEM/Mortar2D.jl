# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

using Base.Test
using Mortar2D: calculate_mortar_matrices, calculate_segments

@testset "calculate bigger interface" begin
    ns = 4 # number of slaves in interface
    nm = 5 # number of masters in interface

    Xs = Dict(i => [x, 1/2] for (i, x) in enumerate(linspace(0, 1, ns)))
    Es = Dict(i => [i, i+1] for i=1:(ns-1))
    Xm = Dict(ns+i => [x, 1/2] for (i, x) in enumerate(linspace(0, 1, nm)))
    Em = Dict(i => [i, i+1]+1 for i=ns:ns+nm-2)

    println(Xs)
    println(Xm)
    println(Es)
    println(Em)

    slave_ids = keys(Es)
    master_ids = keys(Em)
    elements = merge(Es, Em)
    coords = merge(Xs, Xm)
    element_types = Dict(i => :Seg2 for i=1:(ns+nm))
    normals = Dict(i => [0.0, 1.0] for (i, x) in enumerate(linspace(0, 1, ns)))
    segmentation = calculate_segments(slave_ids, master_ids,
                               elements, element_types,
                               coords, normals)
    B = zeros(ns, ns+nm)
    for sid in keys(Es)
        mids = []
        for (mid, xi) in segmentation[sid]
            push!(mids, mid)
        end
        De, Me = calculate_mortar_matrices(sid, mids, elements, element_types,
                                           coords, normals, segmentation)
        B[elements[sid], elements[sid]] += De
        for mid in mids
            B[elements[sid], elements[mid]] += Me[mid]
        end
    end

    # next, calculate projection P
    D = B[1:ns,1:ns]
    M = B[1:ns,ns+1:ns+nm]
    println(D)
    println(M)
    P = inv(D) * M
    P[abs.(P) .< 1.0e-9] = 0
    println(P)
    um = linspace(0, 1, nm)
    us = P*um
    println(us)
    @test isapprox(us, linspace(0, 1, ns))
end

