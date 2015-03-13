This tool is designed to convert triangular (in 2D) or tetrahedral (in 3D) mesh
to quadrilateral or hexahedral one respectively.

At the moment _tethex_ takes a Gmsh's mesh file as input, then converts all 
simplices (tetrahedra and triangles) into bricks (hexahedra and quadrangles), 
and writes the resulting mesh in Gmsh's format again. New mesh can be used 
in many software packages working with bricks only - for example, deal.II - 
the package which this tool was originally designed for.

Supported formats of input meshes:
  * Gmsh ASCII 2.0 (.msh)

Format of output meshes is Gmsh ASCII 2.0 (.msh).

