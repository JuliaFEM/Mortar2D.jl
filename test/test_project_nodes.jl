# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

using Base.Test
using Mortar2D: project_from_slave_to_master, project_from_master_to_slave

@testset "project slave node to master surface" begin
    xm1 = [7.0, 2.0]
    xm2 = [4.0, -2.0]
    xs = [0.0, 0.0]
    ns = [3.0/5.0, -4.0/5.0]
    xi2 = project_from_slave_to_master(Val{:Seg2}, xs, ns, xm1, xm2)
    @test isapprox(xi2, 11/6)
end

@testset "project slave node to master surface" begin
    xm = [4.0, -2.0]
    xs1 = [0.0, 0.0]
    xs2 = [4.0, 3.0]
    ns1 = [3.0/5.0, -4.0/5.0]
    ns2 = sqrt(2)/2*[1, -1]
    xi1 = project_from_master_to_slave(Val{:Seg2}, xm, xs1, xs2, ns1, ns2)
    @test isapprox(xi1, -0.281575016087237)
end
