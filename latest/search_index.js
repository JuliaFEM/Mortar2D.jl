var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Mortar2D.jl-documentation-1",
    "page": "Introduction",
    "title": "Mortar2D.jl documentation",
    "category": "section",
    "text": "(Image: Typical 2d segmentation)Pages = [\"index.md\", \"api.md\"]Mortar2D.jl is a julia package to calculate discrete projections between non-conforming finite element mesheds. The resulting \"mortar matrices\" can be used to tie non-conforming finite element meshes together which are meshed separately to construct bigger models.Using mortar methods in mesh tie problems results variationally consistent solution. Mathematically, goal is to solve mixed problem with primary field variable and Lagrange multipliers, which have a physical meaning (e.g. contact pressure if unknown field is displacement). The problem arising is a typical saddle point problem with zeros on diagonal.Mortar2D.jl is part of JuliaFEM. All codes are MIT licensed."
},

{
    "location": "index.html#Installing-and-testing-package-1",
    "page": "Introduction",
    "title": "Installing and testing package",
    "category": "section",
    "text": "Installing package goes same way like other packages in julia, i.e.julia> Pkg.add(\"Mortar2D\")Testing package can be done using Pkg.test, i.e.julia> Pkg.test(\"Mortar2D\")"
},

{
    "location": "index.html#Contributing-1",
    "page": "Introduction",
    "title": "Contributing",
    "category": "section",
    "text": "Have a new great idea and want to share it with the open source community? From here and here you can look for coding style. Here is explained how to contribute to open source project, in general."
},

{
    "location": "api.html#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api.html#Mortar2D.calculate_normals",
    "page": "API",
    "title": "Mortar2D.calculate_normals",
    "category": "Function",
    "text": "calculate_normals(elements::Dict{Int, Vector{Int}},\n                  element_types::Dict{Int, Symbol},\n                  X::Dict{Int, Vector{Float64})\n\nGiven elements, element types and node locations, calculate nodal normals by first calculating normal directions for each element and then averaging them in nodes. As a result we get unique normal direction defined to each node.\n\nExample\n\nX = Dict(1 => [7.0, 7.0], 2 => [4.0, 3.0], 3 => [0.0, 0.0])\nelements = Dict(1 => [1, 2], 2 => [2, 3])\nelement_types = Dict(1 => :Seg2, 2 => :Seg2)\nnormals = calculate_normals(elements, element_types, X)\n\n# output\n\nDict{Int64,Array{Float64,1}} with 3 entries:\n  2 => [0.707107, -0.707107]\n  3 => [0.6, -0.8]\n  1 => [0.8, -0.6]\n\n\n\n\n"
},

{
    "location": "api.html#Mortar2D.project_from_master_to_slave",
    "page": "API",
    "title": "Mortar2D.project_from_master_to_slave",
    "category": "Function",
    "text": "Find the projection of a master node xm, to the slave surface with nodes (xs1, xs2), in direction of slave surface normal defined by (ns1, ns2).\n\n\n\n"
},

{
    "location": "api.html#Mortar2D.project_from_slave_to_master",
    "page": "API",
    "title": "Mortar2D.project_from_slave_to_master",
    "category": "Function",
    "text": "Find the projection of a slave node xs, having normal vector ns, onto master elements with nodes (xm1, xm2).\n\nMathematics\n\nReferences\n\nYang2005\n\n\n\n"
},

{
    "location": "api.html#Mortar2D.calculate_segments",
    "page": "API",
    "title": "Mortar2D.calculate_segments",
    "category": "Function",
    "text": "Given slave surface elements, master surface elements, nodal coordinates and normal direction on nodes of slave surface elements, calculate contact segments.\n\nReturns a dictionary, where key is slave element id and values is a list of master elements giving contribution to that slave elements and xi-coordinates of slave side element.\n\n\n\n"
},

{
    "location": "api.html#Mortar2D.calculate_mortar_matrices",
    "page": "API",
    "title": "Mortar2D.calculate_mortar_matrices",
    "category": "Function",
    "text": "Given segmentation, calculate mortar matrices De and Me.\n\n\n\n"
},

{
    "location": "api.html#Mortar2D.calculate_mortar_assembly",
    "page": "API",
    "title": "Mortar2D.calculate_mortar_assembly",
    "category": "Function",
    "text": "Given data, calculate projection matrix P. \n\n\n\n"
},

{
    "location": "api.html#API-documentation-1",
    "page": "API",
    "title": "API documentation",
    "category": "section",
    "text": "DocTestSetup = quote\n    using Mortar2D\n    using Mortar2D: calculate_normals, project_from_master_to_slave, project_from_slave_to_master, calculate_segments, calculate_mortar_matrices, calculate_mortar_assembly\nendMortar2D.calculate_normals\nMortar2D.project_from_master_to_slave\nMortar2D.project_from_slave_to_master\nMortar2D.calculate_segments\nMortar2D.calculate_mortar_matrices\nMortar2D.calculate_mortar_assembly"
},

{
    "location": "api.html#Index-1",
    "page": "API",
    "title": "Index",
    "category": "section",
    "text": ""
},

]}
