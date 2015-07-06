/*
 * tethex - tetrahedra to hexahedra conversion
 * Copyright (c) 2013 Mikhail Artemiev
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#include "tethex.h"
#include "config.h"
#include <algorithm>
#include <fstream>
#include <cmath>

TETHEX_NAMESPACE_OPEN

//-------------------------------------------------------
//
// d2s - convert data to string
//
//-------------------------------------------------------
inline std::string d2s(double x, bool scientific, int precision)
{
  std::ostringstream o;
  if (scientific)
  {
    o.setf(std::ios::scientific);
    o.precision(precision);
  }
  if (!(o << x))
    throw std::runtime_error("Bad conversion from double to string!");
  return o.str();
}

inline std::string d2s(int x)
{
  std::ostringstream o;
  if (!(o << x))
    throw std::runtime_error("Bad conversion from int to string!");
  return o.str();
}

inline std::string d2s(unsigned int x)
{
  std::ostringstream o;
  if (!(o << x))
    throw std::runtime_error("Bad conversion from unsigned int to string!");
  return o.str();
}

#if defined(HAVE_64BIT_SIZE_T)
inline std::string d2s(size_t x)
{
  std::ostringstream o;
  if (!(o << x))
    throw std::runtime_error("Bad conversion from size_t to string!");
  return o.str();
}
#endif




//-------------------------------------------------------
//
// expect and require
//
//-------------------------------------------------------
void requirement_fails(const char *file,
                       unsigned int line,
                       std::string message)
{
  std::string exc = "Exception:\nfile = " + std::string(file) +
                    "\nline = " + d2s(line) +
                    "\nmessage = " + message + "\n";
  throw std::runtime_error(exc);
}




//-------------------------------------------------------
//
// Point
//
//-------------------------------------------------------
Point::Point()
{
  for (int i = 0; i < n_coord; ++i)
    coord[i] = 0.;
}

Point::Point(const double coordinates[])
{
  for (int i = 0; i < n_coord; ++i)
    coord[i] = coordinates[i];
}

Point::Point(const double x_coord,
             const double y_coord,
             const double z_coord)
{
  coord[0] = x_coord;
  if (n_coord > 1) coord[1] = y_coord;
  if (n_coord > 2) coord[2] = z_coord;
}

Point::Point(const Point &p)
{
  for (int i = 0; i < n_coord; ++i)
    coord[i] = p.coord[i];
}

Point& Point::operator =(const Point &p)
{
  for (int i = 0; i < n_coord; ++i)
    coord[i] = p.coord[i];
  return *this;
}

double Point::get_coord(unsigned int number) const
{
  expect(number < n_coord,
         "The number of coordinate is incorrect: " +
         d2s(number) + ". It should be in the range: [0, " +
         d2s(n_coord) + ")");

  return coord[number];
}

void Point::set_coord(unsigned int number, double value)
{
  expect(number < n_coord,
         "The number of coordinate is incorrect: " +
         d2s(number) + ". It should be in the range: [0, " +
         d2s(n_coord) + ")");

  coord[number] = value;
}





//-------------------------------------------------------
//
// MeshElement
//
//-------------------------------------------------------
MeshElement::MeshElement(unsigned int n_ver,
                         unsigned int n_edg,
                         unsigned int n_fac,
                         unsigned int el_type)
  : n_vertices(n_ver),
    n_edges(n_edg),
    n_faces(n_fac),
    gmsh_el_type(el_type),
    material_id(0)
{
  vertices.resize(n_vertices, 0);
  edges.resize(n_edges, 0);
  faces.resize(n_faces, 0);
}

MeshElement::~MeshElement()
{
  vertices.clear();
  edges.clear();
  faces.clear();
}

inline unsigned int MeshElement::get_n_vertices() const
{
  expect(n_vertices == vertices.size(),
         "Memory for vertices is not allocated properly (size is " + d2s(vertices.size()) +
         "), or n_vertices (" + d2s(n_vertices) + ") is set to wrong number");
  return n_vertices;
}

inline unsigned int MeshElement::get_n_edges() const
{
  expect(n_edges == edges.size(),
         "Memory for edges is not allocated properly (size is " + d2s(edges.size()) +
         "), or n_edges (" + d2s(n_edges) + ") is set to wrong number");
  return n_edges;
}

inline unsigned int MeshElement::get_n_faces() const
{
  expect(n_faces == faces.size(),
         "Memory for faces is not allocated properly (size is " + d2s(faces.size()) +
         "), or n_faces (" + d2s(n_faces) + ") is set to wrong number");
  return n_faces;
}

inline unsigned int MeshElement::get_gmsh_el_type() const
{
  return gmsh_el_type;
}

inline unsigned int MeshElement::get_material_id() const
{
  return material_id;
}

MeshElement::MeshElement(const MeshElement &elem)
  : n_vertices(elem.n_vertices),
    n_edges(elem.n_edges),
    n_faces(elem.n_faces),
    gmsh_el_type(elem.gmsh_el_type),
    material_id(elem.material_id)
{
  vertices = elem.vertices;
  edges = elem.edges;
  faces = elem.faces;
}

MeshElement& MeshElement::operator =(const MeshElement &elem)
{
  n_vertices = elem.n_vertices;
  n_edges = elem.n_edges;
  n_faces = elem.n_faces;
  gmsh_el_type = elem.gmsh_el_type;
  material_id = elem.material_id;
  vertices = elem.vertices;
  edges = elem.edges;
  faces = elem.faces;
  return *this;
}

unsigned int MeshElement::get_vertex(unsigned int number) const
{
  expect(number < n_vertices,
         "The local number of vertex is incorrect: " + d2s(number) +
         ". It has to be in range [0, " + d2s(n_vertices) + ").");
  return vertices[number];
}

unsigned int MeshElement::get_edge(unsigned int number) const
{
  expect(number < n_edges,
         "The local number of edge is incorrect: " + d2s(number) +
         ". It has to be in range [0, " + d2s(n_edges) + ").");
  return edges[number];
}

unsigned int MeshElement::get_face(unsigned int number) const
{
  expect(number < n_faces,
         "The local number of face is incorrect: " + d2s(number) +
         ". It has to be in range [0, " + d2s(n_faces) + ").");
  return faces[number];
}

void MeshElement::set_vertex(unsigned int local_number, unsigned int global_number)
{
  expect(local_number < get_n_vertices(),
         "Local number (" + d2s(local_number) +
         ") is incorrect. It must be in the range [0, " + d2s(n_edges) + ")");
  vertices[local_number] = global_number;
}

void MeshElement::set_edge(unsigned int local_number, unsigned int global_number)
{
  expect(local_number < get_n_edges(),
         "Local number (" + d2s(local_number) +
         ") is incorrect. It must be in the range [0, " + d2s(n_edges) + ")");
  edges[local_number] = global_number;
}

void MeshElement::set_face(unsigned int local_number, unsigned int global_number)
{
  expect(local_number < get_n_faces(),
         "Local number (" + d2s(local_number) +
         ") is incorrect. It must be in the range [0, " + d2s(n_faces) + ")");
  faces[local_number] = global_number;
}

void MeshElement::set_faces(const std::vector<unsigned int> &face_numbers)
{
  expect(face_numbers.size() == get_n_faces(),
         "Array of face numbers has another size (" + d2s(face_numbers.size()) +
         ") than it must be (" + d2s(get_n_faces()) + ")");
  faces = face_numbers;
}

bool MeshElement::contains(const unsigned int vertex) const
{
  for (unsigned int i = 0; i < n_vertices; ++i)
    if (vertex == vertices[i])
      return true;
  return false;
}




//-------------------------------------------------------
//
// PhysPoint
//
//-------------------------------------------------------
PhysPoint::PhysPoint()
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{ }

PhysPoint::PhysPoint(const std::vector<unsigned int> &ver,
                     const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  expect(vertices.size() == ver.size(),
         "The size of list of vertices (" + d2s(ver.size()) +
         ") is not equal to really needed number of vertices (" + d2s(n_vertices) + ")");
  vertices = ver;
  material_id = mat_id;
}

PhysPoint::PhysPoint(const unsigned int ver,
                     const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  vertices[0] = ver;
  material_id = mat_id;
}




//-------------------------------------------------------
//
// Line
//
//-------------------------------------------------------
Line::Line()
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{ }

Line::Line(const std::vector<unsigned int> &ver,
           const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  expect(vertices.size() == ver.size(),
         "The size of list of vertices (" + d2s(ver.size()) +
         ") is not equal to really needed number of vertices (" + d2s(n_vertices) + ")");
  vertices = ver;
  material_id = mat_id;
}

Line::Line(const unsigned int v1,
           const unsigned int v2,
           const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  vertices[0] = v1;
  vertices[1] = v2;
  material_id = mat_id;
}

unsigned int Line::common_vertex(const Line& line) const
{
  for (unsigned int i = 0; i < n_vertices; ++i)
    if (line.contains(vertices[i]))
      return vertices[i];
  require(false, "There is no common vertex between these two lines!");
  return 0; // to calm compiler down
}

unsigned int Line::another_vertex(const unsigned int vertex) const
{
  if (vertex == vertices[0])
    return vertices[1];
  else if (vertex == vertices[1])
    return vertices[0];
  else
    require(false, "This line doesn't contain the vertex. So we can't find another one.");
  return 0; // to calm compiler down
}





//-------------------------------------------------------
//
// Triangle
//
//-------------------------------------------------------
Triangle::Triangle()
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{ }


Triangle::Triangle(const std::vector<unsigned int> &ver,
                   const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  expect(vertices.size() == ver.size(),
         "The size of list of vertices (" + d2s(ver.size()) +
         ") is not equal to really needed number of vertices (" + d2s(n_vertices) + ")");
  vertices = ver;
  material_id = mat_id;
}

Triangle::Triangle(const unsigned int v1,
                   const unsigned int v2,
                   const unsigned int v3,
                   const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  vertices[0] = v1;
  vertices[1] = v2;
  vertices[2] = v3;
  material_id = mat_id;
}




//-------------------------------------------------------
//
// Tetrahedron
//
//-------------------------------------------------------
Tetrahedron::Tetrahedron()
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{ }


Tetrahedron::Tetrahedron(const std::vector<unsigned int> &ver,
                         const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  expect(vertices.size() == ver.size(),
         "The size of list of vertices (" + d2s(ver.size()) +
         ") is not equal to really needed number of vertices (" + d2s(n_vertices) + ")");
  vertices = ver;
  material_id = mat_id;
}

Tetrahedron::Tetrahedron(const unsigned int v1,
                         const unsigned int v2,
                         const unsigned int v3,
                         const unsigned int v4,
                         const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  vertices[0] = v1;
  vertices[1] = v2;
  vertices[2] = v3;
  vertices[3] = v4;
  material_id = mat_id;
}




//-------------------------------------------------------
//
// Quadrangle
//
//-------------------------------------------------------
Quadrangle::Quadrangle()
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{ }


Quadrangle::Quadrangle(const std::vector<unsigned int> &ver,
                       const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  expect(vertices.size() == ver.size(),
         "The size of list of vertices (" + d2s(ver.size()) +
         ") is not equal to really needed number of vertices (" + d2s(n_vertices) + ")");
  vertices = ver;
  material_id = mat_id;
}

Quadrangle::Quadrangle(const unsigned int v1,
                       const unsigned int v2,
                       const unsigned int v3,
                       const unsigned int v4,
                       const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  vertices[0] = v1;
  vertices[1] = v2;
  vertices[2] = v3;
  vertices[3] = v4;
  material_id = mat_id;
}




//-------------------------------------------------------
//
// Hexahedron
//
//-------------------------------------------------------
Hexahedron::Hexahedron()
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{ }


Hexahedron::Hexahedron(const std::vector<unsigned int> &ver,
                       const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  expect(vertices.size() == ver.size(),
         "The size of list of vertices (" + d2s(ver.size()) +
         ") is not equal to really needed number of vertices (" + d2s(n_vertices) + ")");
  vertices = ver;
  material_id = mat_id;
}

Hexahedron::Hexahedron(const unsigned int v1,
                       const unsigned int v2,
                       const unsigned int v3,
                       const unsigned int v4,
                       const unsigned int v5,
                       const unsigned int v6,
                       const unsigned int v7,
                       const unsigned int v8,
                       const unsigned int mat_id)
  : MeshElement(n_vertices, n_edges, n_faces, gmsh_el_type)
{
  vertices[0] = v1;
  vertices[1] = v2;
  vertices[2] = v3;
  vertices[3] = v4;
  vertices[4] = v5;
  vertices[5] = v6;
  vertices[6] = v7;
  vertices[7] = v8;
  material_id = mat_id;
}





//-------------------------------------------------------
//
// IncidenceMatrix
//
//-------------------------------------------------------
IncidenceMatrix::IncidenceMatrix(const unsigned int n_vertices,
                                 const std::vector<MeshElement_ptr> &cells)
  : dim(n_vertices)
{
  std::vector<unsigned int> *vec = new std::vector<unsigned int>[dim]; // for lower triangle
  // look through all mesh cells
  for (unsigned int cell = 0; cell < cells.size(); ++cell)
  {
    // look at all pairs of cell vertices
    for (unsigned int i = 0; i < cells[cell]->get_n_vertices(); ++i)
    {
      const unsigned int ii = cells[cell]->get_vertex(i);
      for (unsigned int j = 0; j < cells[cell]->get_n_vertices(); ++j)
      {
        const unsigned int jj = cells[cell]->get_vertex(j);
        if (ii > jj) // we consider only lower triangle of matrix
        {
          // add to vector, if the vector doesn't contain this number
          if (std::find(vec[ii].begin(), vec[ii].end(), jj) == vec[ii].end())
            vec[ii].push_back(jj);
        }
      }
    }
  }

  // sorting the vectors
  for (unsigned int i = 0; i < dim; ++i)
    std::sort(vec[i].begin(), vec[i].end());

  // the number of non zero elements in each row of lower triangle
  row = new unsigned int[dim + 1];
  row[0] = 0;
  for (unsigned int i = 0; i < dim; ++i)
    row[i + 1] = row[i] + vec[i].size();

  n_non_zero = row[dim]; // the number of all non zero elements in lower triangle

  // numbers of non zero elements in lower triangle
  col = new unsigned int[n_non_zero];
  unsigned int k = 0;
  for (unsigned int i = 0; i < dim; ++i)
  {
    for (unsigned int j = 0; j < vec[i].size(); ++j)
    {
      col[k] = vec[i][j];
      k++;
    }
  }

  // free the memory
  for (unsigned int i = 0; i < dim; ++i)
    vec[i].clear();
  delete[] vec;
}



IncidenceMatrix::~IncidenceMatrix()
{
  delete[] row;
  delete[] col;
}



unsigned int IncidenceMatrix::find(const unsigned int row_number,
                                   const unsigned int col_number) const
{
  // because we seek in lower triangle, row must be bigger than col
  expect(row_number > col_number,
         "We seek values in lower triangle, so row should be bigger than column. But in this case row_number = " +
         d2s(row_number) + ", col_number = " + d2s(col_number) + "!");

  for (unsigned int i = row[row_number]; i < row[row_number + 1]; ++i)
  {
    if (col[i] == col_number)
      return i;
  }

  // if the number cannot be found in this row, we throw exception
  require(false,
          "The right value wasn't found for such row and column numbers: row_number = " +
          d2s(row_number) + ", col_number = " + d2s(col_number) + "!");
  return 0; // to calm compiler down about returned value
}


inline unsigned int IncidenceMatrix::get_n_nonzero() const
{
  return n_non_zero;
}




//-------------------------------------------------------
//
// Mesh
//
//-------------------------------------------------------
Mesh::Mesh()
  : n_converted_quadrangles(0),
    n_converted_hexahedra(0)
{ }


Mesh::~Mesh()
{
  clean();
}



void Mesh::clean()
{
  vertices.clear();
  points.clear();
  lines.clear();
  edges.clear();
  triangles.clear();
  faces.clear();
  tetrahedra.clear();
  quadrangles.clear();
  hexahedra.clear();
}




void Mesh::read(const std::string &file)
{
  std::ifstream in(file.c_str());
  require(in, "File " + file + " cannot be opened!");

  clean(); // free the memory for mesh elements

  std::string str;
  in >> str; // the first string of Gmsh file is "$MeshFormat"
  expect(str == "$MeshFormat",
         "The first string of the Gmsh file " + file +
         " doesn't equal to \"$MeshFormat\". The actual string is \"" + str + "\"");

  // read the information about the mesh
  double version;
  int binary, dsize;
  in >> version >> binary >> dsize;
  // The function has been testing for meshes corresponding
  // to Gmsh versions since 2.6.0, therefore
  // in debug mode you'll have an exception if the mesh version is less than 2.2.
  // There is no exception in release mode though.
  // So to read the mesh with version 2.1 and less, check that DEBUG variable is set to 0.
  // But NOTE that msh files of 1.0 format have absolutely different structure!
  expect(version >= 2.2,
         "The version of Gmsh's mesh is too old (" + d2s(version) +\
         "). The library was tested for versions 2.2+.");
  expect(dsize == sizeof(double),
         "The size of Gmsh's double (" + d2s(dsize) +\
         ") doesn't equal to size of double type (" + d2s(sizeof(double)) + ")");

  getline(in, str); // read some empty string

  // there is additional 1 (the number - one) in binary format
  if (binary) {
    int one;
    in.read(reinterpret_cast<char*>(&one), sizeof(int));
    require(one == 1,
            "The binary one (" + d2s(one) + ") doesn't equal to 1!");
  }

  // we make a map between serial number of the vertex and its number in the file.
  // it will help us when we create mesh elements
  std::map<unsigned int, unsigned int> vertices_map;

  // read lines of mesh file.
  // if we face specific keyword, we'll treat the section.
  while (in >> str)
  {
    if (str == "$PhysicalNames") // read the section of names of physical entities
    {
      unsigned int n_names;
      in >> n_names;
      getline(in, str);
      physical_names.resize(n_names);
      for (unsigned int i = 0; i < n_names; ++i)
        getline(in, physical_names[i]);
    }

    else if (str == "$Nodes") // read the mesh vertices
    {
      unsigned int n_vertices; // the number of all mesh vertices (that are saved in the file)
      in >> n_vertices; // read that number
      vertices.resize(n_vertices); // allocate the memory for mesh vertices
      getline(in, str); // read some empty string

      unsigned int number; // the number of the vertex
      double coord[Point::n_coord]; // Cartesian coordinates of the vertex (Gmsh produces 3D mesh regardless its real dimension)

      // read vertices
      for (unsigned int ver = 0; ver < n_vertices; ++ver)
      {
        if (binary) // binary format
        {
          in.read(reinterpret_cast<char*>(&number), sizeof(unsigned int)); // the number of each node
          in.read(reinterpret_cast<char*>(coord), Point::n_coord * sizeof(double)); // node coordinates
        }
        else // ASCII format
        {
          in >> number;
          for (unsigned int i = 0; i < Point::n_coord; ++i)
            in >> coord[i];
        }
        vertices[ver] = Point(coord); // save the vertex
        vertices_map[number] = ver; // add the number of vertex to the map
      }

      expect(n_vertices == vertices_map.size(),
             "Vertices numbers are not unique: n_vertices = " + d2s(n_vertices) +
             " vertices_map.size() = " + d2s(vertices_map.size()));

    } // read the vertices

    else if (str == "$Elements") // read the mesh elements
    {
      unsigned int n_elements; // the number of mesh elements
      in >> n_elements; // read that number
      getline(in, str); // empty string

      unsigned int number; // the number of the element [1, nElements]
      unsigned int el_type; // the type of the element (1 - line, 2 - triangle, etc)
      unsigned int n_tags; // the number of tags describing the element
      unsigned int phys_domain; // the physical domain where the element takes place
      unsigned int elem_domain; // the elementary domain where the element takes place
      unsigned int partition; // the partition in which the element takes place

      // the map between the type of the element,
      // and the number of nodes that describe it
      std::map<unsigned int, unsigned int> type_nodes;
      type_nodes[1] = 2; // 2-nodes line
      type_nodes[2] = 3; // 3-nodes triangle
      type_nodes[3] = 4; // 4-nodes quadrangle
      type_nodes[4] = 4; // 4-nodes tetrahedron
      type_nodes[5] = 8; // 8-nodes hexahedron
      type_nodes[15]= 1; // 1-node point

      if (binary) // binary format
      {
/*
        unsigned int n_elem_part = 0; // some part of the elements
        unsigned int header[3]; // the header of the element
        unsigned int amount; // amount of mesh elements of the same type

        while (n_elem_part < n_elements)
        {
          in.read(reinterpret_cast<char*>(header), 3 * sizeof(unsigned int)); // read the header
          el_type = header[0];
          amount  = header[1];
          n_tags  = header[2];

          n_elem_part += amount;

           // the number of nodes
          const unsigned int n_elem_nodes = type_nodes.find(el_type) != ;

          switch (el_type)
          {
          case 1: // 2-node line

            const unsigned int n_data = 1 + n_tags + n_elem_nodes; // how much data we need to read
            int data[n_data]; // allocate memory for data
            unsigned int nodes[n_elem_nodes]; // allocate memory for nodes
            for (unsigned int el = 0; el < amount; ++el)
            {
              in.read(reinterpret_cast<char*>(data), n_data * sizeof(int)); // read the data
              number = data[0];
              phys_domain = (n_tags > 0) ? data[1] : 0; // physical domain - the most important value
              elem_domain = (n_tags > 1) ? data[2] : 0; // elementary domain
              partition   = (n_tags > 2) ? data[3] : 0; // partition (Metis, Chaco, etc)
              for (unsigned int i = 0; i < n_elem_nodes; ++i)
                nodes[i] = vertices_map[data[n_tags + 1 + i]]; // nodes can be numerated not sequentially

              // add new element in the list
              lines.push_back(Edge(nodes[0], nodes[1], phys_domain));
            }
            break;
          case 2: // 3-node triangle
            const unsigned int n_elem_nodes = 3; // the number of nodes
            const unsigned int n_data = 1 + n_tags + n_elem_nodes; // how much data we need to read
            int data[n_data]; // allocate memory for data
            unsigned int nodes[n_elem_nodes]; // allocate memory for nodes
            for (unsigned int el = 0; el < amount; ++el)
            {
              in.read(reinterpret_cast<char*>(data), n_data * sizeof(int)); // read the data
              number = data[0];
              phys_domain = (n_tags > 0) ? data[1] : 0; // physical domain - the most important value
              elem_domain = (n_tags > 1) ? data[2] : 0; // elementary domain
              partition   = (n_tags > 2) ? data[3] : 0; // partition (Metis, Chaco, etc)
              for (unsigned int i = 0; i < n_elem_nodes; ++i)
                nodes[i] = vertices_map[data[n_tags + 1 + i]]; // nodes can be numerated not sequentially

              // add new element in the list
              triangles.push_back(Triangle(nodes[0], nodes[1], nodes[2], phys_domain));
            }
            break;
          case 4: // 4-node tetrahedron
            break;
          default: // other elements are not interesting for us
            break;
          }

          int nElemNodes = mev->second; // the number of nodes
          int nData = 1 + nTags + nElemNodes; // how much data we need to read
          std::vector<int> data(nData); // allocate memory for data
          std::vector<int> nodes(nElemNodes); // allocate memory for nodes

          // read the elements of the same type
          for (int el = 0; el < amount; ++el)
          {
            in.read(reinterpret_cast<char*>(&data[0]), nData * sizeof(int)); // read the data
            number = data[0];
            physDomain = (nTags > 0) ? data[1] : 0; // physical domain - the most important value
            elemDomain = (nTags > 1) ? data[2] : 0; // elementary domain
            partition = (nTags > 2) ? data[3] : 0; // partition (Metis, Chaco, etc)
            for (int i = 0; i < nElemNodes; ++i)
              nodes[i] = data[nTags + 1 + i] - minNodeNumber; // Gmsh numerates the nodes from 'minNodeNumber' (it's usually 1)

            // add new element in the list
            _elements[elName].push_back(MeshElement(nodes, physDomain));
          }

          nodes.clear();
          data.clear();
        }

        // check some expectations
        expect(nElemPart == nElements,
               ExcMessage("The accumulate number of different elements (" + d2s(nElemPart) +\
                          ") is not equal to the amount of all elements in the mesh (" + d2s(nElements) + ")!"));
*/
      } // binary format

      else // ASCII format
      {
        for (int el = 0; el < n_elements; ++el)
        {
          in >> number >> el_type >> n_tags;
          std::vector<unsigned int> data(n_tags); // allocate the memory for some data
          for (unsigned int i = 0; i < n_tags; ++i) // read this information
            in >> data[i];
          phys_domain = (n_tags > 0) ? data[0] : 0; // physical domain - the most important value
          elem_domain = (n_tags > 1) ? data[1] : 0; // elementary domain
          partition   = (n_tags > 2) ? data[2] : 0; // partition (Metis, Chaco, etc)
          data.clear(); // other data isn't interesting for us

          // how many vertices (nodes) describe the element
          std::map<unsigned int, unsigned int>::const_iterator el_type_iter =
              type_nodes.find(el_type);

          require(el_type_iter != type_nodes.end(),
                  "This type of the Gmsh's element (" + d2s(el_type) +
                  ") in the mesh file \"" + file + "\" is unknown!");

          const unsigned int n_elem_nodes = el_type_iter->second; // the number of nodes
          std::vector<unsigned int> nodes(n_elem_nodes); // allocate memory for nodes
          for (unsigned int i = 0; i < n_elem_nodes; ++i)
          {
            in >> nodes[i]; // read the numbers of nodes
            // vertices can be numerated not sequentially (or not from 0)
            nodes[i] = vertices_map.find(nodes[i])->second;
          }

          // add new element in the list
          //MeshElement *new_element;
          switch (el_type)
          {
          case 15: // 1-node point
            points.push_back(MeshElement_ptr(new PhysPoint(nodes, phys_domain)));
            break;
          case 1: // 2-nodes line
            lines.push_back(MeshElement_ptr(new Line(nodes, phys_domain)));
            //new_element = new Line(nodes, phys_domain);
            break;
          case 2: // 3-nodes triangle
            triangles.push_back(MeshElement_ptr(new Triangle(nodes, phys_domain)));
            //new_element = new Triangle(nodes, phys_domain);
            break;
          case 3: // 4-nodes quadrangle
            quadrangles.push_back(MeshElement_ptr(new Quadrangle(nodes, phys_domain)));
            //new_element = new Quadrangle(nodes, phys_domain);
            break;
          case 4: //4-nodes tetrahedron
            tetrahedra.push_back(MeshElement_ptr(new Tetrahedron(nodes, phys_domain)));
            //new_element = new Tetrahedron(nodes, phys_domain);
            break;
          case 5: // 8-nodes hexahedron
            hexahedra.push_back(MeshElement_ptr(new Hexahedron(nodes, phys_domain)));
            //new_element = new Hexahedron(nodes, phys_domain);
            break;
          default:
            require(false,
                    "Unknown type of the Gmsh's element (" + d2s(el_type) +
                    ") in the file " + file + "!");
          }

          nodes.clear();

          //elements.push_back(new_element);
          //delete new_element;
        }

        // check some expectations
        expect(number == n_elements,
               "The number of the last read Gmsh's element (" + d2s(number) +\
               ") is not equal to the amount of all elements in the mesh (" + d2s(n_elements) + ")!");

      } // ASCII format

      // requirements after reading elements
      require(!triangles.empty() || !tetrahedra.empty() ||
              !quadrangles.empty() || !hexahedra.empty(),
             "There are no any 2D or 3D elements in the mesh!");

      // we prevent mixing of simplices and bricks in one mesh.
      // at least for the moment
      if (!triangles.empty() || !tetrahedra.empty())
        require(quadrangles.empty() && hexahedra.empty(),
                "There are simplices and bricks in the same mesh. It's prohibited. Mesh file " + file);

    } // read the elements
  }

  in.close(); // close the file
}




void Mesh::convert()
{
  if (!tetrahedra.empty())
    convert_3D();
  else if (!triangles.empty())
    convert_2D();

  if (!hexahedra.empty())
    convert_hexahedra();

  if (!quadrangles.empty())
    convert_quadrangles();
}



void Mesh::set_new_vertices(const std::vector<MeshElement_ptr> &elements,
                            const unsigned int n_old_vertices,
                            const unsigned int shift)
{
  for (unsigned int elem = 0; elem < elements.size(); ++elem)
  {
    for (unsigned int coord = 0; coord < Point::n_coord; ++coord)
    {
      double coordinate = 0.;
      for (unsigned int ver = 0; ver < elements[elem]->get_n_vertices(); ++ver)
      {
        const unsigned int cur_vertex = elements[elem]->get_vertex(ver);
        expect(cur_vertex < n_old_vertices,
               "The element has a vertex (" + d2s(cur_vertex) +
               ") that is more than the number of old vertices (" + d2s(n_old_vertices) + ")");
        coordinate += vertices[cur_vertex].get_coord(coord);
      }
      vertices[n_old_vertices + shift + elem].set_coord(coord,
                                                        coordinate / elements[elem]->get_n_vertices());
    }
  }
}




void Mesh::convert_2D()
{
  // firstly we need to numerate all edges
  const IncidenceMatrix incidence_matrix(vertices.size(), triangles);

  // third parameter - whether we need to initialize the vector of all edges of the mesh.
  // yes - we need it
  edge_numeration(triangles, incidence_matrix, true);

  // after edge numbering
  // we should add new nodes -
  // one node at the middle of every edge and
  // one node at the center of every triangle
  const unsigned int n_old_vertices = vertices.size();
  vertices.resize(n_old_vertices + edges.size() + triangles.size());

  // add 'edge'-nodes - at the middle of edge
  set_new_vertices(edges, n_old_vertices, 0);
  //add_edge_nodes(n_old_vertices);

  // add 'triangle'-nodes - at the center of triangle
  set_new_vertices(triangles, n_old_vertices, edges.size());
  //add_triangle_nodes(triangles, n_old_vertices);

  // convert triangles into quadrangles.
  // third parameter - whether we need to numerate edges of triangles,
  // no - we've already done it
  convert_triangles(incidence_matrix, n_old_vertices, false);

  // now we don't need triangles anymore
#if defined(DELETE_SIMPLICES)
  triangles.clear();
#endif

  // after that we check boundary elements (lines),
  // because after adding new vertices they need to be redefined
  redefine_lines(incidence_matrix, n_old_vertices);

#if !defined(TESTING)
  edges.clear();
#endif
}




void Mesh::convert_3D()
{
  // firstly - edge numeration
  const IncidenceMatrix incidence_matrix(vertices.size(), tetrahedra);

  // third parameter - whether we need to initialize the vector of all edges of the mesh.
  // yes - we need it
  edge_numeration(tetrahedra, incidence_matrix, true);

  // the main structure - this vector of maps
  // vector's element is associated with mesh edges and vector
  // has size = n_edges.
  // the map's element - pair of numbers:
  // the key - number of vertex, opposite to the suitable edge,
  // the value - number of suitable face, defining by edge and vertex
  // opposite to edge
  VectorMap edge_vertex_incidence(edges.size());

  // after that - face numeration
  face_numeration(tetrahedra, incidence_matrix, edge_vertex_incidence);

  // some checks
  expect(vertices.size() + faces.size() - 1 == tetrahedra.size() + edges.size(),
         "Some sophisticated assumption is not held!");

  // after edge and face numbering
  // we should add new nodes -
  // one node at the middle of every edge,
  // one node at the center of every face and
  // one node at the center of every tetrahedron
  const unsigned int n_old_vertices = vertices.size();
  vertices.resize(n_old_vertices + edges.size() + faces.size() + tetrahedra.size());

  // add 'edge'-nodes - at the middle of edge
  set_new_vertices(edges, n_old_vertices, 0);
  //add_edge_nodes(n_old_vertices);

  // add 'face'-nodes - at the center of face
  set_new_vertices(faces, n_old_vertices, edges.size());
  //add_triangle_nodes(faces, n_old_vertices);

  // add 'tetrahedron'-nodes - at the center of every tetrahedron
  set_new_vertices(tetrahedra, n_old_vertices, edges.size() + faces.size());

  // now we generate hexahedrons
  convert_tetrahedra(n_old_vertices, incidence_matrix, edge_vertex_incidence);

  // now we don't need tetrahedra anymore
#if defined(DELETE_SIMPLICES)
  tetrahedra.clear();
#endif

  // third parameter - whether we need to numerate edges of triangles,
  // yes, because for boundary triangles we haven't done it before
  convert_triangles(incidence_matrix, n_old_vertices, true, edge_vertex_incidence);

  // now we don't need triangles anymore
#if defined(DELETE_SIMPLICES)
  triangles.clear();
#endif

  // after that we check lines (1D boundary elements),
  // because after adding new vertices they need to be redefined
  redefine_lines(incidence_matrix, n_old_vertices);

#if !defined(TESTING)
  edges.clear();
  faces.clear();
#endif

}




void Mesh::convert_tetrahedra(const unsigned int n_old_vertices,
                              const IncidenceMatrix &incidence_matrix,
                              const VectorMap edge_vertex_incidence)
{
  std::vector<unsigned int> hexahedron_vertices(Hexahedron::n_vertices);

  for (unsigned int tet = 0; tet < tetrahedra.size(); ++tet)
  {
    for (unsigned int ver = 0; ver < Tetrahedron::n_vertices; ++ver)
    {
      // current vertex
      const unsigned int cur_vertex = tetrahedra[tet]->get_vertex(ver);

      // we're looking for 3 edges to which this vertex belongs
      std::vector<unsigned int> seek_edges;
      for (unsigned int edge = 0; edge < Tetrahedron::n_edges; ++edge)
      {
        const unsigned int cur_edge = tetrahedra[tet]->get_edge(edge);
        if (edges[cur_edge]->contains(cur_vertex))
          seek_edges.push_back(cur_edge);
      }
      expect(seek_edges.size() == 3, "");

      // numeration of hexahedron vertices
      hexahedron_vertices[0] = cur_vertex;
      hexahedron_vertices[1] = n_old_vertices + seek_edges[0];
      hexahedron_vertices[2] = n_old_vertices + edges.size() +
                               find_face_from_two_edges(seek_edges[0], seek_edges[1],
                                                        incidence_matrix, edge_vertex_incidence);
      hexahedron_vertices[3] = n_old_vertices + seek_edges[1];
      hexahedron_vertices[4] = n_old_vertices + seek_edges[2];
      hexahedron_vertices[5] = n_old_vertices + edges.size() +
                               find_face_from_two_edges(seek_edges[0], seek_edges[2],
                                                        incidence_matrix, edge_vertex_incidence);
      hexahedron_vertices[6] = n_old_vertices + edges.size() + faces.size() + tet;
      hexahedron_vertices[7] = n_old_vertices + edges.size() +
                               find_face_from_two_edges(seek_edges[1], seek_edges[2],
                                                        incidence_matrix, edge_vertex_incidence);

      seek_edges.clear();

      // check cell measure to be sure that we numerate all hexahedra in one way.
      // this measure taken from deal.II.
      change_vertices_order(3, vertices, hexahedron_vertices);

//      // convert the order of vertices to suitable for deal.II to check the cell measure
//      std::vector<unsigned int> vertices_dealII_order(Hexahedron::n_vertices);
//      unsigned int order_to_deal[] = { 0, 1, 5, 4, 2, 3, 7, 6 };

//      for (unsigned int i = 0; i < Hexahedron::n_vertices; ++i)
//        vertices_dealII_order[order_to_deal[i]] = hexahedron_vertices[i];

//      if (cell_measure_3D(vertices, vertices_dealII_order) < 0)
//        // reorder vertices - swap front and back faces
//        for (unsigned int i = 0; i < Quadrangle::n_vertices; ++i)
//          std::swap(hexahedron_vertices[i], hexahedron_vertices[i + 4]);

      // now generate hexahedron
      hexahedra.push_back(MeshElement_ptr(new Hexahedron(hexahedron_vertices,
                                                         tetrahedra[tet]->get_material_id())));

    } // vertices
  } // tetrahedra

  require(tetrahedra.size() * 4 == hexahedra.size(),
          "The number of hexahedra (" + d2s(hexahedra.size()) +
          ") is not equal to number of tetrahedra (" + d2s(tetrahedra.size()) +
          ") multiplying by 3 (" + d2s(3 * tetrahedra.size()) + ")");
}




void Mesh::convert_triangles(const IncidenceMatrix &incidence_matrix,
                             const unsigned int n_old_vertices,
                             bool numerate_edges,
                             const VectorMap &edge_vertex_incidence)
{
  // quadrangles generation
  std::vector<unsigned int> quadrangle_vertices(Quadrangle::n_vertices);

  // we need to numerate edges of boundary triangles
  // numerate_edges = true - the case of 3D mesh, when we numerate edges of boundary triangles,
  //                         but in this case we mustn't numerate edges themselves,
  //                         because they were already numerated during process of numeration edges of tetrahedra,
  // numerate_edges = false - the case of 2D mesh, when we have already numerated edges of triangles,
  //                          so we don't need to do that again
  if (numerate_edges)
    // third parameter - whether we need to initialize the vector of all edges of the mesh,
    // no we shouldn't
    edge_numeration(triangles, incidence_matrix, false);

  for (unsigned int tri = 0; tri < triangles.size(); ++tri)
  {
    for (unsigned int ver = 0; ver < Triangle::n_vertices; ++ver)
    {
      // current vertex
      const unsigned int cur_vertex = triangles[tri]->get_vertex(ver);

      // we are looking for 2 edges that contain current vertex
      std::vector<unsigned int> seek_edges;
      for (unsigned int edge = 0; edge < Triangle::n_edges; ++edge)
      {
        const unsigned int cur_edge = triangles[tri]->get_edge(edge);
        if (edges[cur_edge]->contains(cur_vertex))
          seek_edges.push_back(cur_edge);

      }
      expect(seek_edges.size() == 2,
             "The number of edges to which every vertex belongs must be equal to 2");

      // numeration of quadrangle vertices
      quadrangle_vertices[0] = cur_vertex;
      quadrangle_vertices[1] = n_old_vertices + seek_edges[0];
      // !!! need to repair !!!
      if (numerate_edges) // distinguish 2D and 3D cases
      {
        // 3D case - boundary triangles
        quadrangle_vertices[2] = n_old_vertices + edges.size() +
                                 find_face_from_two_edges(seek_edges[0],
                                                          seek_edges[1],
                                                          incidence_matrix,
                                                          edge_vertex_incidence);
      }
      else
      {
        // 2D case - triangles themselves
        quadrangle_vertices[2] = n_old_vertices + edges.size() + tri;
      }
      quadrangle_vertices[3] = n_old_vertices + seek_edges[1];

      seek_edges.clear();

      // though the order of vertices is right, it may be clockwise or counterclockwise,
      // and it's important not to mix these 2 directions.
      // so, we need additional check as deal.II authors do.
      change_vertices_order(2, vertices, quadrangle_vertices);

//      // taken from deal.II
//      if (cell_measure_2D(vertices, quadrangle_vertices) < 0)
//        // change 2 vertices to reverse the order
//        std::swap(quadrangle_vertices[1], quadrangle_vertices[3]);

      // now we are ready to generate quadrangle
      quadrangles.push_back(MeshElement_ptr(new Quadrangle(quadrangle_vertices,
                                                           triangles[tri]->get_material_id())));

    } // for every vertex we have one quadrangle
  } // triangles

  require(triangles.size() * 3 == quadrangles.size(),
          "The number of quadrangles (" + d2s(quadrangles.size()) +
          ") is not equal to number of triangles (" + d2s(triangles.size()) +
          ") multiplying by 3 (" + d2s(3 * triangles.size()) + ")");
}




void Mesh::convert_quadrangles()
{
  for (unsigned int elem = 0; elem < quadrangles.size(); ++elem)
  {
    std::vector<unsigned int> quad_vertices(Quadrangle::n_vertices);
    for (unsigned int i = 0; i < Quadrangle::n_vertices; ++i)
      quad_vertices[i] = quadrangles[elem]->get_vertex(i);

    const unsigned ver = quad_vertices[1];

    change_vertices_order(2, vertices, quad_vertices);

    // since only first and third vertices are swapped (if any)
    // we compare only one vertex number
    if (quad_vertices[1] != ver)
      ++n_converted_quadrangles;

    for (unsigned int i = 0; i < Quadrangle::n_vertices; ++i)
      quadrangles[elem]->set_vertex(i, quad_vertices[i]);
  }
}





void Mesh::convert_hexahedra()
{
  for (unsigned int elem = 0; elem < hexahedra.size(); ++elem)
  {
    std::vector<unsigned int> hexahedron_vertices(Hexahedron::n_vertices);
    for (unsigned int i = 0; i < Hexahedron::n_vertices; ++i)
      hexahedron_vertices[i] = hexahedra[elem]->get_vertex(i);

    const unsigned ver = hexahedron_vertices[0];

    change_vertices_order(3, vertices, hexahedron_vertices);

    // we compare only one vertex number
    if (hexahedron_vertices[0] != ver)
      ++n_converted_hexahedra;

    for (unsigned int i = 0; i < Hexahedron::n_vertices; ++i)
      hexahedra[elem]->set_vertex(i, hexahedron_vertices[i]);
  }
}




void Mesh::redefine_lines(const IncidenceMatrix &incidence_matrix,
                          const unsigned int n_old_vertices)
{
  const unsigned int n_old_lines = lines.size();
  for (unsigned int line = 0; line < n_old_lines; ++line)
  {
    // we need to find an edge that coincides with this line
    const unsigned int ver1 = lines[line]->get_vertex(0);
    const unsigned int ver2 = lines[line]->get_vertex(1);
    const unsigned int edge = incidence_matrix.find(std::max(ver1, ver2),
                                                    std::min(ver1, ver2));

    // we change existing line and add new line at the end of list
    lines[line]->set_vertex(1, n_old_vertices + edge); // changing existing line
    lines.push_back(MeshElement_ptr(new Line(n_old_vertices + edge, ver2,
                                             lines[line]->get_material_id()))); // add new line
  }

  require(n_old_lines * 2 == lines.size(),
          "The number of physical lines (" + d2s(lines.size()) +
          ") is not equal to number of original physical lines (" + d2s(n_old_lines) +
          ") multiplying by 2 (" + d2s(2 * n_old_lines) + ")");
}





void Mesh::edge_numeration(std::vector<MeshElement_ptr> &cells,
                           const IncidenceMatrix &incidence_matrix,
                           bool initialize_edges)
{
  //// matrix of incidence between vertices of the mesh
  //const IncidenceMatrix incidence_matrix(vertices.size(), cells);

  // the number of edges in such mesh - the number of non zero elements in incidence matrix
  const unsigned int n_edges = incidence_matrix.get_n_nonzero();

  // allocate memory for edges
  if (initialize_edges)
    edges.resize(n_edges);

  // look through all cells of the mesh
  for (unsigned int cell = 0; cell < cells.size(); ++cell)
  {
    unsigned int lne = 0; // local number of the edge (0 <= lne < cell::n_edges)
    for (unsigned int i = 0; i < cells[cell]->get_n_vertices(); ++i)
    {
      const unsigned int ii = cells[cell]->get_vertex(i);
      for (unsigned int j = 0; j < cells[cell]->get_n_vertices(); ++j)
      {
        const unsigned int jj = cells[cell]->get_vertex(j);
        if (ii > jj) // ii must be bigger than jj
        {
          const unsigned int gne = incidence_matrix.find(ii, jj); // global number of edge
          // set the global number of edge to cell
          cells[cell]->set_edge(lne, gne);
          // initialize edge
          if (initialize_edges)
            edges[gne] = MeshElement_ptr(new Line(std::min(ii, jj),
                                                  std::max(ii, jj),
                                                  cells[cell]->get_material_id()));
          // increase local number of edge
          ++lne;
        }
      }
    }
    expect(lne == cells[cell]->get_n_edges(),
           "lne must be equal to " + d2s(cells[cell]->get_n_edges()) +
           ", but it is " + d2s(lne));
  }

} // edge numeration




void Mesh::face_numeration(std::vector<MeshElement_ptr> &cells,
                           const IncidenceMatrix &incidence_matrix,
                           VectorMap &edge_vertex_incidence)
{
  unsigned int n_faces = 0; // the number of all faces and the number of current face

  for (unsigned int cell = 0; cell < cells.size(); ++cell)
  {
    std::vector<unsigned int> face_numbers;
    for (unsigned int edge = 0; edge < cells[cell]->get_n_edges(); ++edge)
    {
      unsigned int cur_edge = cells[cell]->get_edge(edge);
      for (unsigned int ver = 0; ver < cells[cell]->get_n_vertices(); ++ver)
      {
        unsigned int cur_vertex = cells[cell]->get_vertex(ver);
        if (!edges[cur_edge]->contains(cur_vertex)) // if edge doesn't contain vertex - they are opposite to each other
        {
          // edge and vertex opposite to it - they define a face
          if (edge_vertex_incidence[cur_edge].find(cur_vertex) == edge_vertex_incidence[cur_edge].end())
          {
            // if there is no such pair of edge and vertex -
            // that means that this face was not numerated yet.
            // so do it now
            edge_vertex_incidence[cur_edge][cur_vertex] = n_faces;

            // and we should do it for all pairs of edges and opposite vertices
            // for this face to avoid duplicates.
            unsigned int another_edge = incidence_matrix.find(std::max(edges[cur_edge]->get_vertex(0), cur_vertex),
                                                              std::min(edges[cur_edge]->get_vertex(0), cur_vertex));
            edge_vertex_incidence[another_edge][edges[cur_edge]->get_vertex(1)] = n_faces;

            // and once more time
            another_edge = incidence_matrix.find(std::max(edges[cur_edge]->get_vertex(1), cur_vertex),
                                                 std::min(edges[cur_edge]->get_vertex(1), cur_vertex));
            edge_vertex_incidence[another_edge][edges[cur_edge]->get_vertex(0)] = n_faces;

            // create array of face's vertices
            std::vector<unsigned int> face_vertices(3);
            face_vertices[0] = edges[cur_edge]->get_vertex(0);
            face_vertices[1] = edges[cur_edge]->get_vertex(1);
            face_vertices[2] = cur_vertex;
            sort(face_vertices.begin(), face_vertices.end());

            // add this face to faces list
            faces.push_back(MeshElement_ptr(new Triangle(face_vertices,
                                                         cells[cell]->get_material_id())));

            // add the number of face into array
            face_numbers.push_back(n_faces);

            // increase the number of faces
            ++n_faces;
          }
          else // this face already exsists
          {
            // add the number of face into array
            // if there is no such number there yet
            const unsigned int fnumber = edge_vertex_incidence[cur_edge].find(cur_vertex)->second;
            if (find(face_numbers.begin(), face_numbers.end(), fnumber) == face_numbers.end())
              face_numbers.push_back(fnumber);
          }
        } // find opposite vertex
      } // vertices
    } // edges

    expect(face_numbers.size() == cells[cell]->get_n_faces(),
           "There is no enough faces for " + d2s(cell) +
           "-th cell. It's " + d2s(face_numbers.size()) +
           ". But is must be " + d2s(cells[cell]->get_n_faces()));

    // set these face numbers as cell's faces
    cells[cell]->set_faces(face_numbers);

  } // cells
} // face numeration





void Mesh::write(const std::string &file)
{
  std::ofstream out(file.c_str());
  require(out, "File " + file + " cannot be opened for writing!");

  out.setf(std::ios::scientific);
  out.precision(16);

  out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";
  if (!physical_names.empty())
  {
    out << "$PhysicalNames\n";
    out << physical_names.size() << "\n";
    for (unsigned int i = 0; i < physical_names.size(); ++i)
      out << physical_names[i] << "\n";
    out << "$EndPhysicalNames\n";
  }

  out << "$Nodes\n" << vertices.size() << "\n";
  for (unsigned int ver = 0; ver < vertices.size(); ++ver)
  {
    out << ver + 1 << " ";
    for (unsigned int coord = 0; coord < Point::n_coord; ++coord)
      out << vertices[ver].get_coord(coord) << " ";
    out << "\n";
  }

  const unsigned int n_all_elements = points.size() +
                                      lines.size() +
                                      triangles.size() +
                                      tetrahedra.size() +
                                      quadrangles.size() +
                                      hexahedra.size();
  out << "$EndNodes\n$Elements\n" << n_all_elements << "\n";

  unsigned int serial_number = 0;

  write_elements(out, points, serial_number);
  write_elements(out, lines, serial_number);
  write_elements(out, triangles, serial_number);
  write_elements(out, tetrahedra, serial_number);
  write_elements(out, quadrangles, serial_number);
  write_elements(out, hexahedra, serial_number);

  out << "$EndElements\n";

  out.close();
}



unsigned int Mesh::get_n_vertices() const
{
  return vertices.size();
}

unsigned int Mesh::get_n_points() const
{
  return points.size();
}

unsigned int Mesh::get_n_lines() const
{
  return lines.size();
}

unsigned int Mesh::get_n_edges() const
{
  return edges.size();
}

unsigned int Mesh::get_n_triangles() const
{
  return triangles.size();
}

unsigned int Mesh::get_n_faces() const
{
  return faces.size();
}

unsigned int Mesh::get_n_tetrahedra() const
{
  return tetrahedra.size();
}

unsigned int Mesh::get_n_quadrangles() const
{
  return quadrangles.size();
}

unsigned int Mesh::get_n_hexahedra() const
{
  return hexahedra.size();
}



Point Mesh::get_vertex(const unsigned int number) const
{
  expect(number < vertices.size(),
         "The required number (" + d2s(number) +
         " is bigger than the number of vertices (" + d2s(vertices.size()) + "))");
  return vertices[number];
}

MeshElement& Mesh::get_point(const unsigned int number) const
{
  expect(number < points.size(),
         "The required number (" + d2s(number) +
         " is bigger than the number of physical points (" + d2s(points.size()) + "))");
  return *(points[number]);
}

MeshElement& Mesh::get_edge(const unsigned int number) const
{
  expect(number < edges.size(),
         "The required number (" + d2s(number) +
         " is bigger than the number of edges (" + d2s(edges.size()) + "))");
  return *(edges[number]);
}

MeshElement& Mesh::get_line(const unsigned int number) const
{
  expect(number < lines.size(),
         "The required number (" + d2s(number) +
         " is bigger than the number of lines (" + d2s(lines.size()) + "))");
  return *(lines[number]);
}

MeshElement& Mesh::get_face(const unsigned int number) const
{
  expect(number < faces.size(),
         "The required number (" + d2s(number) +
         " is bigger than the number of faces (" + d2s(faces.size()) + "))");
  return *(faces[number]);
}

MeshElement& Mesh::get_triangle(const unsigned int number) const
{
  expect(number < triangles.size(),
         "The required number (" + d2s(number) +
         " is bigger than the number of triangles (" + d2s(triangles.size()) + "))");
  return *(triangles[number]);
}

MeshElement& Mesh::get_tetrahedron(const unsigned int number) const
{
  expect(number < tetrahedra.size(),
         "The required number (" + d2s(number) +
         " is bigger than the number of tetrahedra (" + d2s(tetrahedra.size()) + "))");
  return *(tetrahedra[number]);
}

MeshElement& Mesh::get_quadrangle(const unsigned int number) const
{
  expect(number < quadrangles.size(),
         "The required number (" + d2s(number) +
         " is bigger than the number of quadrangles (" + d2s(quadrangles.size()) + "))");
  return *(quadrangles[number]);
}

MeshElement& Mesh::get_hexahedron(const unsigned int number) const
{
  expect(number < hexahedra.size(),
         "The required number (" + d2s(number) +
         " is bigger than the number of hexahedra (" + d2s(hexahedra.size()) + "))");
  return *(hexahedra[number]);
}



void Mesh::info(std::ostream &out) const
{
  out << "\nvertices       : " << vertices.size()
      << "\npoints (phys)  : " << points.size()
      << "\nedges          : " << edges.size()
      << "\nlines          : " << lines.size()
      << "\ntriangles      : " << triangles.size()
      << "\nfaces          : " << faces.size()
      << "\ntetrahedra     : " << tetrahedra.size()
      << "\nquadrangles    : " << quadrangles.size()
      << "\nhexahedra      : " << hexahedra.size()
      << "\nconverted quads: " << n_converted_quadrangles
      << "\nconverted hexs : " << n_converted_hexahedra
      << "\n\n";
}



void Mesh::statistics(std::ostream &out) const
{
  out << "\nvertices:\n";
  for (unsigned int i = 0; i < vertices.size(); ++i)
  {
    out << i << " ";
    for (unsigned int j = 0; j < Point::n_coord; ++j)
      out << vertices[i].get_coord(j) << " ";
    out << "\n";
  }

  out << "\nedges:\n";
  for (unsigned int i = 0; i < edges.size(); ++i)
  {
    out << i << " ";
    for (unsigned int j = 0; j < Line::n_vertices; ++j)
      out << edges[i]->get_vertex(j) << " ";
    out << "\n";
  }

  out << "\ntets:\n";
  for (unsigned int i = 0; i < tetrahedra.size(); ++i)
  {
    out << i << "  ";
    for (unsigned int j = 0; j < Tetrahedron::n_vertices; ++j)
      out << tetrahedra[i]->get_vertex(j) << " ";
    out << " | ";
    for (unsigned int j = 0; j < Tetrahedron::n_edges; ++j)
      out << tetrahedra[i]->get_edge(j) << " ";
    out << " | ";
    for (unsigned int j = 0; j < Tetrahedron::n_faces; ++j)
      out << tetrahedra[i]->get_face(j) << " ";
    out << "\n";
  }

  out << "\nfaces:\n";
  for (unsigned int i = 0; i < faces.size(); ++i)
  {
    out << i << " ";
    for (unsigned int j = 0; j < faces[i]->get_n_vertices(); ++j)
      out << faces[i]->get_vertex(j) << " ";
    out << "\n";
  }
  out << "\n";

}



unsigned int Mesh::find_face_from_two_edges(const unsigned int edge1,
                                            const unsigned int edge2,
                                            const IncidenceMatrix &vertices_incidence,
                                            const VectorMap &edge_vertex_incidence) const
{
  // initialize auxiliary lines
  Line line1(edges[edge1]->get_vertex(0), edges[edge1]->get_vertex(1), edges[edge1]->get_material_id());
  Line line2(edges[edge2]->get_vertex(0), edges[edge2]->get_vertex(1), edges[edge2]->get_material_id());

  // find common vertex
  const unsigned int common_vertex = line1.common_vertex(line2);

  // find other 2 vertices
  const unsigned int ver1 = line1.another_vertex(common_vertex);
  const unsigned int ver2 = line2.another_vertex(common_vertex);

  // find opposite edge
  const unsigned int opposite_edge = vertices_incidence.find(std::max(ver1, ver2), std::min(ver1, ver2));

  // find the number of face
  return edge_vertex_incidence[opposite_edge].find(common_vertex)->second;
}





//-------------------------------------------------------
//
// Auxiliary functions
//
//-------------------------------------------------------
void write_elements(std::ostream &out,
                    const std::vector<MeshElement_ptr> &elems,
                    unsigned int &serial_number)
{
  for (unsigned int el = 0; el < elems.size(); ++el, ++serial_number)
  {
    out << serial_number + 1 << " "             /* serial number of element */
        << elems[el]->get_gmsh_el_type()        /* type of element suitable for Gmsh */
        << " 2 "                                /* the number of tags */
        << elems[el]->get_material_id() << " "  /* physical domain */
        << elems[el]->get_material_id() << " "; /* elemetary domain - let it be the same */
    for (unsigned int ver = 0; ver < elems[el]->get_n_vertices(); ++ver)
      out << elems[el]->get_vertex(ver) + 1 << " ";
    out << "\n";
  }
}




void change_vertices_order(int dimension,
                           const std::vector<Point> &all_mesh_vertices,
                           std::vector<unsigned int> &vertices)
{
  if (dimension == 2)
  {
    // convert the order of vertices to suitable for deal.II to check the cell measure
    std::vector<unsigned int> vertices_dealII_order(Quadrangle::n_vertices);
    unsigned int order_to_deal[] = { 0, 1, 3, 2 };

    for (unsigned int i = 0; i < Quadrangle::n_vertices; ++i)
      vertices_dealII_order[order_to_deal[i]] = vertices[i];

    if (cell_measure_2D(all_mesh_vertices, vertices_dealII_order) < 0)
      // reorder vertices - swap first and third vertices
      std::swap(vertices[1], vertices[3]);
  }
  else if (dimension == 3)
  {
    // convert the order of vertices to suitable for deal.II to check the cell measure
    std::vector<unsigned int> vertices_dealII_order(Hexahedron::n_vertices);
    unsigned int order_to_deal[] = { 0, 1, 5, 4, 2, 3, 7, 6 };

    for (unsigned int i = 0; i < Hexahedron::n_vertices; ++i)
      vertices_dealII_order[order_to_deal[i]] = vertices[i];

    if (cell_measure_3D(all_mesh_vertices, vertices_dealII_order) < 0)
      // reorder vertices - swap front and back faces
      for (unsigned int i = 0; i < 4; ++i)
        std::swap(vertices[i], vertices[i + 4]);
  }
  else
    require(false, "This feature is not implemented!");
}




double cell_measure_2D(const std::vector<Point> &vertices,
                       const std::vector<unsigned int> &indices)
{
  const double x[] = { vertices[indices[0]].get_coord(0),
                       vertices[indices[1]].get_coord(0),
                       vertices[indices[2]].get_coord(0),
                       vertices[indices[3]].get_coord(0)
                     };
  const double y[] = { vertices[indices[0]].get_coord(1),
                       vertices[indices[1]].get_coord(1),
                       vertices[indices[2]].get_coord(1),
                       vertices[indices[3]].get_coord(1)
                     };
  return (-x[1]*y[0]+x[1]*y[3]+
          y[0]*x[2]+x[0]*y[1]-
          x[0]*y[2]-y[1]*x[3]-
          x[2]*y[3]+x[3]*y[2]) / 2.;
}



double cell_measure_3D(const std::vector<Point> &vertices,
                       const std::vector<unsigned int> &indices)
{
  const double x[8] = { vertices[indices[0]].get_coord(0),
                        vertices[indices[1]].get_coord(0),
                        vertices[indices[2]].get_coord(0),
                        vertices[indices[3]].get_coord(0),
                        vertices[indices[4]].get_coord(0),
                        vertices[indices[5]].get_coord(0),
                        vertices[indices[6]].get_coord(0),
                        vertices[indices[7]].get_coord(0)
                      };
  const double y[8] = { vertices[indices[0]].get_coord(1),
                        vertices[indices[1]].get_coord(1),
                        vertices[indices[2]].get_coord(1),
                        vertices[indices[3]].get_coord(1),
                        vertices[indices[4]].get_coord(1),
                        vertices[indices[5]].get_coord(1),
                        vertices[indices[6]].get_coord(1),
                        vertices[indices[7]].get_coord(1)
                      };
  const double z[8] = { vertices[indices[0]].get_coord(2),
                        vertices[indices[1]].get_coord(2),
                        vertices[indices[2]].get_coord(2),
                        vertices[indices[3]].get_coord(2),
                        vertices[indices[4]].get_coord(2),
                        vertices[indices[5]].get_coord(2),
                        vertices[indices[6]].get_coord(2),
                        vertices[indices[7]].get_coord(2)
                      };

  const double t3 = y[3]*x[2];
  const double t5 = z[1]*x[5];
  const double t9 = z[3]*x[2];
  const double t11 = x[1]*y[0];
  const double t14 = x[4]*y[0];
  const double t18 = x[5]*y[7];
  const double t20 = y[1]*x[3];
  const double t22 = y[5]*x[4];
  const double t26 = z[7]*x[6];
  const double t28 = x[0]*y[4];
  const double t34 = z[3]*x[1]*y[2]+t3*z[1]-t5*y[7]+y[7]*x[4]*z[6]+t9*y[6]-
                     t11*z[4]-t5*y[3]-t14*z[2]+z[1]*x[4]*y[0]-t18*z[3]+
                     t20*z[0]-t22*z[0]-y[0]*x[5]*z[4]-t26*y[3]+t28*z[2]-
                     t9*y[1]-y[1]*x[4]*z[0]-t11*z[5];
  const double t37 = y[1]*x[0];
  const double t44 = x[1]*y[5];
  const double t46 = z[1]*x[0];
  const double t49 = x[0]*y[2];
  const double t52 = y[5]*x[7];
  const double t54 = x[3]*y[7];
  const double t56 = x[2]*z[0];
  const double t58 = x[3]*y[2];
  const double t64 = -x[6]*y[4]*z[2]-t37*z[2]+t18*z[6]-x[3]*y[6]*z[2]+
                     t11*z[2]+t5*y[0]+t44*z[4]-t46*y[4]-t20*z[7]-t49*z[6]-
                     t22*z[1]+t52*z[3]-t54*z[2]-t56*y[4]-t58*z[0]+
                     y[1]*x[2]*z[0]+t9*y[7]+t37*z[4];
  const double t66 = x[1]*y[7];
  const double t68 = y[0]*x[6];
  const double t70 = x[7]*y[6];
  const double t73 = z[5]*x[4];
  const double t76 = x[6]*y[7];
  const double t90 = x[4]*z[0];
  const double t92 = x[1]*y[3];
  const double t95 = -t66*z[3]-t68*z[2]-t70*z[2]+t26*y[5]-t73*y[6]-
                     t14*z[6]+t76*z[2]-t3*z[6]+x[6]*y[2]*z[4]-
                     z[3]*x[6]*y[2]+t26*y[4]-t44*z[3]-x[1]*y[2]*z[0]+
                     x[5]*y[6]*z[4]+t54*z[5]+t90*y[2]-t92*z[2]+t46*y[2];
  const double t102 = x[2]*y[0];
  const double t107 = y[3]*x[7];
  const double t114 = x[0]*y[6];
  const double t125 = y[0]*x[3]*z[2]-z[7]*x[5]*y[6]-x[2]*y[6]*z[4]+
                      t102*z[6]-t52*z[6]+x[2]*y[4]*z[6]-t107*z[5]-
                      t54*z[6]+t58*z[6]-x[7]*y[4]*z[6]+t37*z[5]-t114*z[4]+
                      t102*z[4]-z[1]*x[2]*y[0]+t28*z[6]-y[5]*x[6]*z[4]-z[5]*x[1]*y[4]-t73*y[7];
  const double t129 = z[0]*x[6];
  const double t133 = y[1]*x[7];
  const double t145 = y[1]*x[5];
  const double t156 = t90*y[6]-t129*y[4]+z[7]*x[2]*y[6]-t133*z[5]+
                      x[5]*y[3]*z[7]-t26*y[2]-t70*z[3]+t46*y[3]+
                      z[5]*x[7]*y[4]+z[7]*x[3]*y[6]-t49*z[4]+t145*z[7]-
                      x[2]*y[7]*z[6]+t70*z[5]+t66*z[5]-z[7]*x[4]*y[6]+t18*z[4]+x[1]*y[4]*z[0];
  const double t160 = x[5]*y[4];
  const double t165 = z[1]*x[7];
  const double t178 = z[1]*x[3];
  const double t181 = t107*z[6]+t22*z[7]+t76*z[3]+t160*z[1]-x[4]*y[2]*z[6]+
                      t70*z[4]+t165*y[5]+x[7]*y[2]*z[6]-t76*z[5]-t76*z[4]+
                      t133*z[3]-t58*z[1]+y[5]*x[0]*z[4]+t114*z[2]-
                      t3*z[7]+t20*z[2]+t178*y[7]+t129*y[2];
  const double t207 = t92*z[7]+t22*z[6]+z[3]*x[0]*y[2]-x[0]*y[3]*z[2]-
                      z[3]*x[7]*y[2]-t165*y[3]-t9*y[0]+t58*z[7]+
                      y[3]*x[6]*z[2]+t107*z[2]+t73*y[0]-x[3]*y[5]*z[7]+
                      t3*z[0]-t56*y[6]-z[5]*x[0]*y[4]+t73*y[1]-t160*z[6]+t160*z[0];
  const double t228 = -t44*z[7]+z[5]*x[6]*y[4]-t52*z[4]-t145*z[4]+t68*z[4]+
                      t92*z[5]-t92*z[0]+t11*z[3]+t44*z[0]+t178*y[5]-t46*y[5]-
                      t178*y[0]-t145*z[0]-t20*z[5]-t37*z[3]-
                      t160*z[7]+t145*z[3]+x[4]*y[6]*z[2];

  return (t34+t64+t95+t125+t156+t181+t207+t228)/12.;
}


TETHEX_NAMESPACE_CLOSE
