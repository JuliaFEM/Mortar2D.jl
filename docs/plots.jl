# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

using PyPlot
using Mortar2D: calculate_normals, project_from_master_to_slave, project_from_slave_to_master, calculate_segments, calculate_mortar_matrices, calculate_mortar_assembly

coords = Dict(1 => [8.0, 10.0],
         2 => [7.0, 7.0],
         3 => [4.0, 3.0],
         4 => [0.0, 0.0],
         5 => [-3.0, 0.0],
         6 => [12.0, 10.0],
         7 => [10.0, 4.0],
         8 => [7.0, 2.0],
         9 => [4.0, -2.0],
         10 => [0.0, -3.0],
         11 => [-4.0, -3.0])

elements = Dict(
    1 => [1, 2],
    2 => [2, 3],
    3 => [3, 4],
    4 => [4, 5],
    5 => [6, 7],
    6 => [7, 8],
    7 => [8, 9],
    8 => [9, 10],
    9 => [10, 11])

slave_element_ids = [1, 2, 3, 4]
slave_elements = Dict(i => elements[i] for i in slave_element_ids)
master_element_ids = [5, 6, 7, 8, 9]
element_types = Dict(i => :Seg2 for i=1:length(elements))
normals = calculate_normals(slave_elements, element_types, coords)
segmentation = calculate_segments(slave_element_ids, master_element_ids,
                                  elements, element_types, coords, normals)

function plot1(; plot_element_normals=false, plot_nodal_normals=false, plot_segmentation=false)

    figure(figsize=(4, 4))

    for i in slave_element_ids
        x1,y1 = coords[elements[i][1]]
        x2,y2 = coords[elements[i][2]]
        xmid = 1/2*(x1+x2)
        ymid = 1/2*(y1+y2)
        n = [y1-y2, x2-x1]
        n /= norm(n)
        plot([x1,x2], [y1,y2], "-bo")
        if plot_element_normals
            arrow(x1, y1, n[1], n[2], head_width=0.1, head_length=0.2, fc="b", ec="b")
            arrow(x2, y2, n[1], n[2], head_width=0.1, head_length=0.2, fc="b", ec="b")
        end
    end

    for i in master_element_ids
        x1,y1 = coords[elements[i][1]]
        x2,y2 = coords[elements[i][2]]
        plot([x1,x2], [y1,y2], "-ro")
    end

    if plot_nodal_normals
        for i in keys(normals)
            x = coords[i]
            n = normals[i]
            arrow(x[1], x[2], n[1], n[2], head_width=0.1, head_length=0.2, fc="k", ec="k")
        end
    end

    if plot_segmentation
        for (sid, seg) in segmentation
            scon = elements[sid]
            xs1 = coords[scon[1]]
            xs2 = coords[scon[2]]
            ns1 = normals[scon[1]]
            ns2 = normals[scon[2]]
            for (mid, xi) in seg
                mcon = elements[mid]
                xm1 = coords[mcon[1]]
                xm2 = coords[mcon[2]]
                for xi_s in xi
                    #xi_s = (1-xi)/2*p + (1+xi)/2*p
                    N1 = [(1-xi_s)/2 (1+xi_s)/2]
                    n_s = N1[1]*ns1 + N1[2]*ns2
                    x_g = N1[1]*xs1 + N1[2]*xs2
                    xi_m = project_from_slave_to_master(Val{:Seg2}, x_g, n_s, xm1, xm2)
                    N2 = [(1-xi_m)/2 (1+xi_m)/2]
                    x_m = N2[1]*xm1 + N2[2]*xm2
                    plot([x_g[1], x_m[1]], [x_g[2], x_m[2]], "--k")
                end
            end
        end
    end

    axis("equal")
    annotate("\$\\Gamma_{\\mathrm{c}}^{\\left(1\\right)}\$", xy=(2.5, 5.0))
    annotate("\$\\Gamma_{\\mathrm{c}}^{\\left(2\\right)}\$", xy=(8.0, 0.0))
    axis("off")
end
