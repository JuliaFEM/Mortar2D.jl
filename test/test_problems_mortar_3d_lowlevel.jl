# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

using JuliaFEM
using JuliaFEM.Preprocess
using JuliaFEM.Postprocess
using JuliaFEM.Testing

@testset "forget to add elements to problem" begin
    X = Dict(
        1 => [0.0, 0.0, 0.0],
        2 => [1.0, 0.0, 0.0],
        3 => [0.0, 1.0, 0.0],
        4 => [-0.25, 0.50, 0.00],
        5 => [0.50, -0.25, 0.00],
        6 => [0.75, 0.75, 0.00])
    s = Element(Tri3, [1, 2, 3])
    m = Element(Tri3, [4, 5, 6])
    update!([s, m], "geometry", X)
    update!(s, "master elements", [m])
    p = Problem(Mortar, "two elements", 1, "temperature")
    initialize!(p)
    assemble!(p)
    @test true
end

""" Calculate mortar projection matrix P = D^-1*M from mortar assembly. """
function calculate_mortar_projection_matrix(problem::Problem{Mortar}, ndim::Int)

    C1 = sparse(problem.assembly.C1, ndim, ndim)
    C2 = sparse(problem.assembly.C2, ndim, ndim)

    @assert nnz(sparse(problem.assembly.K)) == 0
    @assert nnz(sparse(problem.assembly.D)) == 0
    @assert nnz(sparse(problem.assembly.Kg)) == 0
    @assert nnz(sparse(problem.assembly.fg)) == 0
    @assert nnz(sparse(problem.assembly.f)) == 0
    @assert nnz(sparse(problem.assembly.g)) == 0

    @assert C1 == C2
    #@assert problem.properties.dual_basis == true
    @assert problem.properties.adjust == false
 
    S = get_nonzero_rows(C2)
    M = setdiff(get_nonzero_columns(C2), S)

    # Construct matrix P = D^-1*M
    D_ = C2[S,S]
    M_ = -C2[S,M]

    #=
    P = nothing
    if !isdiag(D_)
        warn("D is not diagonal, is dual basis used? This might take a long time.")
        P = ldltfact(1/2*(D_ + D_')) \ M_
    else
        P = D_ \ M_
    end
    =#

    P = lufact(D_) \ full(M_)

    return S, M, P
end

@testset "two linear element clipping, calculation of projection matrix P for standard and dual basis" begin
    X = Dict(
        1 => [0.0, 0.0, 0.0],
        2 => [1.0, 0.0, 0.0],
        3 => [0.0, 1.0, 0.0],
        4 => [-0.25, 0.50, 0.00],
        5 => [0.50, -0.25, 0.00],
        6 => [0.75, 0.75, 0.00])
    s = Element(Tri3, [1, 2, 3])
    m = Element(Tri3, [4, 5, 6])
    update!([s, m], "geometry", X)
    update!(s, "master elements", [m])
    p = Problem(Mortar, "two elements", 1, "temperature")
    p.properties.dual_basis = false
    p.elements = [s; m]
    initialize!(p)
    assemble!(p)
    C1 = sparse(p.assembly.C1)
    C2 = sparse(p.assembly.C2)
    D = sparse(p.assembly.D)
    @test length(D) == 0
    @test C1 == C2
    S, M, P = calculate_mortar_projection_matrix(p, 6)
    @test S == [1, 2, 3]
    @test M == [4, 5, 6]
    # visually inspected to be ok result
    P_expected = 1/15*[9 9 -3; -7 13 9; 13 -7 9]
    @test isapprox(P, P_expected)
    um = [7.5, 15.0, 22.5]
    @test isapprox(P*um, [9.0, 23.0, 13.0])

    empty!(p.assembly)
    p.properties.dual_basis = true
    assemble!(p)
    C1 = sparse(p.assembly.C1)
    C2 = sparse(p.assembly.C2)
    D = sparse(p.assembly.D)
    @test length(D) == 0
    @test C1 == C2
    S, M, P = calculate_mortar_projection_matrix(p, 6)
    @test S == [1, 2, 3]
    @test M == [4, 5, 6]
    @test isapprox(P, P_expected)
    um = [7.5, 15.0, 22.5]
    @test isapprox(P*um, [9.0, 23.0, 13.0])
end


@testset "two quadratic element clipping, calculation of projection matrix P for standard basis" begin
    X = Dict(
        1 => [0.0, 0.0, 0.0],
        2 => [1.0, 0.0, 0.0],
        3 => [0.0, 1.0, 0.0],
        7 => [-0.25, 0.50, 0.00],
        8 => [0.50, -0.25, 0.00],
        9 => [0.75, 0.75, 0.00])
    # middle nodes
    X[4] = 1/2*(X[1] + X[2])
    X[5] = 1/2*(X[2] + X[3])
    X[6] = 1/2*(X[3] + X[1])
    X[10] = 1/2*(X[7] + X[8])
    X[11] = 1/2*(X[8] + X[9])
    X[12] = 1/2*(X[9] + X[7])
    s = Element(Tri6, [1, 2, 3, 4, 5, 6])
    m = Element(Tri6, [7, 8, 9, 10, 11, 12])
    update!([s, m], "geometry", X)
    update!(s, "master elements", [m])
    p = Problem(Mortar, "two elements", 1, "temperature")
    p.properties.dual_basis = false
    p.properties.alpha = 0.2
    p.elements = [s; m]
    initialize!(p)
    assemble!(p)
    C1 = sparse(p.assembly.C1)
    C2 = sparse(p.assembly.C2)
    D = sparse(p.assembly.D)
    @test length(D) == 0
    @test C1 == C2
    S, M, P = calculate_mortar_projection_matrix(p, 12)
    @test S == [1, 2, 3, 4, 5, 6]
    @test M == [7, 8, 9, 10, 11, 12]
    println(full(P))
    # visually inspected to be ok result
    P_expected = 1/675*[81 81 189 972 -324 -324; 609 429 81 -1092 1404 -756; 429 609 81 -1092 -756 1404; -39 231 -81 132 396 36; -81 -81 81 108 324 324; 231 -39 -81 132 36 396]
    @test isapprox(P, P_expected)
    um = 15/2*[1, 2, 3]
    um = [um[1], um[2], um[3], 0.5*(um[1]+um[2]), 0.5*(um[2]+um[3]), 0.5*(um[3]+um[1])]
    us = P*um
    @test isapprox(us, [9.0, 23.0, 13.0, 16.0, 18.0, 11.0])
end

@testset "two quadratic element clipping, calculation of projection matrix P for dual lagrange basis" begin
    X = Dict(
        1 => [0.0, 0.0, 0.0],
        2 => [1.0, 0.0, 0.0],
        3 => [0.0, 1.0, 0.0],
        7 => [-0.25, 0.50, 0.00],
        8 => [0.50, -0.25, 0.00],
        9 => [0.75, 0.75, 0.00])
    # middle nodes
    X[4] = 1/2*(X[1] + X[2])
    X[5] = 1/2*(X[2] + X[3])
    X[6] = 1/2*(X[3] + X[1])
    X[10] = 1/2*(X[7] + X[8])
    X[11] = 1/2*(X[8] + X[9])
    X[12] = 1/2*(X[9] + X[7])
    s = Element(Tri6, [1, 2, 3, 4, 5, 6])
    m = Element(Tri6, [7, 8, 9, 10, 11, 12])
    update!([s, m], "geometry", X)
    update!(s, "master elements", [m])
    p = Problem(Mortar, "two elements", 1, "temperature")
    p.properties.dual_basis = true
    p.properties.alpha = 0.2
    p.elements = [s; m]
    initialize!(p)
    assemble!(p)
    C1 = sparse(p.assembly.C1)
    C2 = sparse(p.assembly.C2)
    D = sparse(p.assembly.D)
    @test length(D) == 0
    @test C1 == C2
    S, M, P = calculate_mortar_projection_matrix(p, 12)
    P_expected = 1/675*[81 81 189 972 -324 -324; 609 429 81 -1092 1404 -756; 429 609 81 -1092 -756 1404; -39 231 -81 132 396 36; -81 -81 81 108 324 324; 231 -39 -81 132 36 396]
    @test isapprox(P, P_expected)
    um = 15/2*[1, 2, 3]
    um = [um[1], um[2], um[3], 0.5*(um[1]+um[2]), 0.5*(um[2]+um[3]), 0.5*(um[3]+um[1])]
    us = P*um
    @test isapprox(us, [9.0, 23.0, 13.0, 16.0, 18.0, 11.0])
end

