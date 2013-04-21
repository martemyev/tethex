#include "tethex.h"
#include <algorithm>
#include <fstream>
#include <map>

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




////-------------------------------------------------------
////
//// Edge
////
////-------------------------------------------------------
//Edge::Edge()
//{
//  vertices[0] = vertices[1] = 0;
//}

//Edge::Edge(unsigned int v1, unsigned int v2)
//{
//  vertices[0] = std::min(v1, v2);
//  vertices[1] = std::max(v1, v2);
//}

//Edge::Edge(const Edge& e)
//{
//  vertices[0] = e.vertices[0];
//  vertices[1] = e.vertices[1];
//}

//Edge& Edge::operator =(const Edge& e)
//{
//  vertices[0] = e.vertices[0];
//  vertices[1] = e.vertices[1];
//  return *this;
//}

//inline unsigned int Edge::get_beg() const
//{
//  return vertices[0];
//}

//inline unsigned int Edge::get_end() const
//{
//  return vertices[1];
//}




//-------------------------------------------------------
//
// Line
//
//-------------------------------------------------------
Line::Line()
  : material_id(0)
{
  vertices[0] = vertices[1] = 0;
}


Line::Line(const std::vector<unsigned int> &ver, const unsigned int mat_id)
  : material_id(mat_id)
{
  expect(ver.size() == n_vertices,
         "The size of vector of vertices (" +
         d2s(ver.size()) + ") is not equal to " + d2s(n_vertices) + "!");
  vertices[0] = ver[0];
  vertices[1] = ver[1];
}

Line::Line(unsigned int v1, unsigned int v2, unsigned int mat_id)
  : material_id(mat_id)
{
  vertices[0] = v1;
  vertices[1] = v2;
}

Line::Line(const Line& line)
  : material_id(line.material_id)
{
  vertices[0] = line.vertices[0];
  vertices[1] = line.vertices[1];
}

Line& Line::operator =(const Line& line)
{
  vertices[0] = line.vertices[0];
  vertices[1] = line.vertices[1];
  material_id = line.material_id;
  return *this;
}

//bool Line::operator ==(const Edge &edge) const
//{
//  return ((vertices[0] == edge.get_beg() && vertices[1] == edge.get_end())
//          ||
//          (vertices[0] == edge.get_end() && vertices[1] == edge.get_beg()));
//}

bool Line::operator ==(const Line &line) const
{
  return ((vertices[0] == line.get_beg() && vertices[1] == line.get_end())
          ||
          (vertices[0] == line.get_end() && vertices[1] == line.get_beg()));
}


inline unsigned int Line::get_beg() const
{
  return vertices[0];
}

inline unsigned int Line::get_end() const
{
  return vertices[1];
}

inline unsigned int Line::get_material_id() const
{
  return material_id;
}

inline void Line::set_end(unsigned int ver)
{
  vertices[1] = ver;
}





//-------------------------------------------------------
//
// Triangle
//
//-------------------------------------------------------
Triangle::Triangle()
  : material_id(0)
{
  for (unsigned int i = 0; i < n_vertices; ++i)
    vertices[i] = 0;
  for (unsigned int i = 0; i < n_edges; ++i)
    edges[i] = 0;
}


Triangle::Triangle(const std::vector<unsigned int> &ver, const unsigned int mat_id)
  : material_id(mat_id)
{
  expect(ver.size() == n_vertices,
         "The size of vector of vertices (" +
         d2s(ver.size()) + ") is not equal to " + d2s(n_vertices) + "!");
  for (unsigned int i = 0; i < n_vertices; ++i)
    vertices[i] = ver[i];
}


Triangle::Triangle(const Triangle &tri)
  : material_id(tri.material_id)
{
  for (unsigned int i = 0; i < n_vertices; ++i)
    vertices[i] = tri.vertices[i];
  for (unsigned int i = 0; i < n_edges; ++i)
    edges[i] = tri.edges[i];
}



Triangle& Triangle::operator =(const Triangle &tri)
{
  for (unsigned int i = 0; i < n_vertices; ++i)
    vertices[i] = tri.vertices[i];
  for (unsigned int i = 0; i < n_edges; ++i)
    edges[i] = tri.edges[i];
  material_id = tri.material_id;
  return *this;
}



inline unsigned int Triangle::get_vertex(unsigned int number) const
{
  expect(number < n_vertices,
         "The local number of vertex is incorrect: " + d2s(number) +
         ". It has to be in range [0, " + d2s(n_vertices) + ").");
  return vertices[number];
}



inline unsigned int Triangle::get_edge(unsigned int number) const
{
  expect(number < n_edges,
         "The local number of edge is incorrect: " + d2s(number) +
         ". It has to be in range [0, " + d2s(n_edges) + ").");
  return edges[number];
}



inline unsigned int Triangle::get_material_id() const
{
  return material_id;
}


void Triangle::set_edge(unsigned int local_number, unsigned int global_number)
{
  expect(local_number < n_edges,
         "Local number (" + d2s(local_number) +
         ") is incorrect. It must be in the range [0, " + d2s(n_edges) + ")");
  edges[local_number] = global_number;
}




//-------------------------------------------------------
//
// Tetrahedron
//
//-------------------------------------------------------
Tetrahedron::Tetrahedron()
  : material_id(0)
{
  for (unsigned int i = 0; i < n_vertices; ++i)
    vertices[i] = 0;
  for (unsigned int i = 0; i < n_edges; ++i)
    edges[i] = 0;
}


Tetrahedron::Tetrahedron(const std::vector<unsigned int> &ver, const unsigned int mat_id)
  : material_id(mat_id)
{
  expect(ver.size() == n_vertices,
         "The size of vector of vertices (" +
         d2s(ver.size()) + ") is not equal to " + d2s(n_vertices) + "!");
  for (unsigned int i = 0; i < n_vertices; ++i)
    vertices[i] = ver[i];
}


Tetrahedron::Tetrahedron(const Tetrahedron &tet)
  : material_id(tet.material_id)
{
  for (unsigned int i = 0; i < n_vertices; ++i)
    vertices[i] = tet.vertices[i];
  for (unsigned int i = 0; i < n_edges; ++i)
    edges[i] = tet.edges[i];
}



Tetrahedron& Tetrahedron::operator =(const Tetrahedron &tet)
{
  for (unsigned int i = 0; i < n_vertices; ++i)
    vertices[i] = tet.vertices[i];
  for (unsigned int i = 0; i < n_edges; ++i)
    edges[i] = tet.edges[i];
  material_id = tet.material_id;
  return *this;
}



inline unsigned int Tetrahedron::get_vertex(unsigned int number) const
{
  expect(number < n_vertices,
         "The local number of vertex is incorrect: " + d2s(number) +
         ". It has to be in range [0, " + d2s(n_vertices) + ").");
  return vertices[number];
}



inline unsigned int Tetrahedron::get_edge(unsigned int number) const
{
  expect(number < n_edges,
         "The local number of edge is incorrect: " + d2s(number) +
         ". It has to be in range [0, " + d2s(n_edges) + ").");
  return edges[number];
}



inline unsigned int Tetrahedron::get_material_id() const
{
  return material_id;
}




//-------------------------------------------------------
//
// Quadrangle
//
//-------------------------------------------------------
Quadrangle::Quadrangle()
  : material_id(0)
{
  for (unsigned int i = 0; i < n_vertices; ++i)
    vertices[i] = 0;
  //for (unsigned int i = 0; i < n_edges; ++i)
  //  edges[i] = 0;
}


Quadrangle::Quadrangle(const std::vector<unsigned int> &ver, const unsigned int mat_id)
  : material_id(mat_id)
{
  expect(ver.size() == n_vertices,
         "The size of vector of vertices (" +
         d2s(ver.size()) + ") is not equal to " + d2s(n_vertices) + "!");
  for (unsigned int i = 0; i < n_vertices; ++i)
    vertices[i] = ver[i];
}


Quadrangle::Quadrangle(const Quadrangle &quad)
  : material_id(quad.material_id)
{
  for (unsigned int i = 0; i < n_vertices; ++i)
    vertices[i] = quad.vertices[i];
  //for (unsigned int i = 0; i < n_edges; ++i)
  //  edges[i] = quad.edges[i];
}



Quadrangle& Quadrangle::operator =(const Quadrangle &quad)
{
  for (unsigned int i = 0; i < n_vertices; ++i)
    vertices[i] = quad.vertices[i];
  //for (unsigned int i = 0; i < n_edges; ++i)
  //  edges[i] = quad.edges[i];
  material_id = quad.material_id;
  return *this;
}



inline unsigned int Quadrangle::get_vertex(unsigned int number) const
{
  expect(number < n_vertices,
         "The local number of vertex is incorrect: " + d2s(number) +
         ". It has to be in range [0, " + d2s(n_vertices) + ").");
  return vertices[number];
}



//inline unsigned int Quadrangle::get_edge(unsigned int number) const
//{
//  expect(number < n_edges,
//         "The local number of edge is incorrect: " + d2s(number) +
//         ". It has to be in range [0, " + d2s(n_edges) + ").");
//  return edges[number];
//}



inline unsigned int Quadrangle::get_material_id() const
{
  return material_id;
}


//void Quadrangle::set_edge(unsigned int local_number, unsigned int global_number)
//{
//  expect(local_number < n_edges,
//         "Local number (" + d2s(local_number) +
//         ") is incorrect. It must be in the range [0, " + d2s(n_edges) + ")");
//  edges[local_number] = global_number;
//}






//-------------------------------------------------------
//
// IncidenceMatrix
//
//-------------------------------------------------------
IncidenceMatrix::IncidenceMatrix(const unsigned int n_vertices,
                                 const std::vector<Triangle> &cells)
  : dim(n_vertices)
{
  std::vector<unsigned int> *vec = new std::vector<unsigned int>[dim]; // for lower triangle
  // look through all mesh cells
  for (unsigned int cell = 0; cell < cells.size(); ++cell)
  {
    // look at all pairs of cell vertices
    for (unsigned int i = 0; i < Triangle::n_vertices; ++i)
    {
      const unsigned int ii = cells[cell].get_vertex(i);
      for (unsigned int j = 0; j < Triangle::n_vertices; ++j)
      {
        const unsigned int jj = cells[cell].get_vertex(j);
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
{
}


Mesh::~Mesh()
{
  clean();
}



void Mesh::clean()
{
  vertices.clear();
  lines.clear();
  edges.clear();
  triangles.clear();
  tetrahedra.clear();
  quadrangles.clear();
}




void Mesh::read(const std::string &file)
{
  std::ifstream in(file.c_str());
  require(in, "File " + file + " cannot be opened!");

  clean(); // free the memory for mesh elements

  std::string str;
  getline(in, str); // the first string of Gmsh file is "$MeshFormat"
  expect(str == "$MeshFormat",
         "The first string of the Gmsh file " + file +
         " doesn't equal to \"$MeshFormat\". The actual string is " + str);

  // read the information about the mesh
  double version;
  int binary, dsize;
  in >> version >> binary >> dsize;
  // The function has been testing for meshes corresponding
  // to Gmsh versions since 2.6.0, therefore
  // in debug mode you'll have an exception if the mesh version is less than 2.2.
  // There is no exception in release mode though.
  // So to read the mesh with version 2.1-, check that DEBUG variable is set to 0.
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
  while (getline(in, str))
  {
    if (str == "$Nodes") // read the mesh vertices
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
      type_nodes[4] = 4; // 4-nodes tetrahedron
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
                  "This type of the element (" + d2s(el_type) +
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
          switch (el_type)
          {
          case 1: // 2-nodes line
            lines.push_back(Line(nodes, phys_domain));
            break;
          case 2: // 3-nodes triangle
            triangles.push_back(Triangle(nodes, phys_domain));
            break;
          case 4: //4-nodes tetrahedron
            tetrahedra.push_back(Tetrahedron(nodes, phys_domain));
            break;
          default:
            require(false,
                    "Unknown type of the element (" + d2s(el_type) +
                    ") in the file " + file + "!");
          }

          nodes.clear();
        }

        // check some expectations
        expect(number == n_elements,
               "The number of the last read element (" + d2s(number) +\
               ") is not equal to the amount of all elements in the mesh (" + d2s(n_elements) + ")!");
      } // ASCII format

      // expectations after reading elements
      expect(!lines.empty() || !triangles.empty() || !tetrahedra.empty(),
             "There is no elements in the mesh!");

    } // read the elements
  }

  in.close(); // close the file
}




void Mesh::convert()
{
  // firstly we need to numerate all edges
  edge_numeration();

  // after edge numbering
  // we should add new nodes -
  // one node at the middle of every edge and
  // one node at the center of every triangle
  const unsigned int n_old_vertices = vertices.size();
  vertices.resize(n_old_vertices + edges.size() + triangles.size());

  // add 'edge'-nodes - at the middle of edge
  for (unsigned int edge = 0; edge < edges.size(); ++edge)
  {
    const unsigned int beg_ver = edges[edge].get_beg();
    const unsigned int end_ver = edges[edge].get_end();
    expect(beg_ver < n_old_vertices,
           "The first vertex (" + d2s(beg_ver) +
           ") of edge (" + d2s(edge) + ") is more than number of vertices (" +
           d2s(n_old_vertices) + ")");
    expect(end_ver < n_old_vertices,
           "The second vertex (" + d2s(end_ver) +
           ") of edge (" + d2s(edge) + ") is more than number of vertices (" +
           d2s(n_old_vertices) + ")");
    for (unsigned int coord = 0; coord < Point::n_coord; ++coord)
    {
      vertices[n_old_vertices + edge].set_coord(coord,
                                                (vertices[beg_ver].get_coord(coord) +
                                                 vertices[end_ver].get_coord(coord)) / 2.);
    }
  }

  // add 'triangle'-nodes - at the center of triangle
  for (unsigned int tri = 0; tri < triangles.size(); ++tri)
  {
    for (unsigned int coord = 0; coord < Point::n_coord; ++coord)
    {
      double coordinate = 0.;
      for (unsigned int ver = 0; ver < Triangle::n_vertices; ++ver)
      {
        expect(triangles[tri].get_vertex(ver) < n_old_vertices, "");
        coordinate += vertices[triangles[tri].get_vertex(ver)].get_coord(coord);
      }
      vertices[n_old_vertices + edges.size() + tri].set_coord(coord, coordinate / Triangle::n_vertices);
    }
  }

  // now we generate quadrangles
  std::vector<unsigned int> quadrangle_vertices(Quadrangle::n_vertices);

  for (unsigned int tri = 0; tri < triangles.size(); ++tri)
  {
    for (unsigned int ver = 0; ver < Triangle::n_vertices; ++ver)
    {
      const unsigned int vertex_number = triangles[tri].get_vertex(ver);
      quadrangle_vertices[0] = vertex_number;
      int number_of_quad_vert = 1;
      for (unsigned int edge = 0; edge < Triangle::n_edges; ++edge)
      {
        const unsigned int edge_number = triangles[tri].get_edge(edge);
        if ((vertex_number == edges[edge_number].get_beg())
            ||
            (vertex_number == edges[edge_number].get_end()))
          quadrangle_vertices[number_of_quad_vert++] = n_old_vertices + edge_number;
      }
      expect(number_of_quad_vert - 1 == 2,
             "The number of edges to which every vertex belongs must be equal to 2");

      quadrangle_vertices[number_of_quad_vert] = n_old_vertices + edges.size() + tri;
      // swap 2 values to make right order
      std::swap(quadrangle_vertices[2], quadrangle_vertices[3]);

      // but though the order of vertices is right it may be clockwise and counterclockwise,
      // and it's important not to mix these 2 directions.
      // so, we need additional check as deal.II authors do.
      // taken from deal.II
      const double x[] = { vertices[quadrangle_vertices[0]].get_coord(0),
                           vertices[quadrangle_vertices[1]].get_coord(0),
                           vertices[quadrangle_vertices[2]].get_coord(0),
                           vertices[quadrangle_vertices[3]].get_coord(0)
                         };
      const double y[] = { vertices[quadrangle_vertices[0]].get_coord(1),
                           vertices[quadrangle_vertices[1]].get_coord(1),
                           vertices[quadrangle_vertices[2]].get_coord(1),
                           vertices[quadrangle_vertices[3]].get_coord(1)
                         };
      const double cell_measure = (-x[1]*y[0]+x[1]*y[3]+
                                    y[0]*x[2]+x[0]*y[1]-
                                    x[0]*y[2]-y[1]*x[3]-
                                    x[2]*y[3]+x[3]*y[2]) / 2.;
      if (cell_measure < 0)
        std::swap(quadrangle_vertices[1], quadrangle_vertices[3]);

      // now we are ready to generate quadrangle
      Quadrangle quad(quadrangle_vertices, triangles[tri].get_material_id());
      quadrangles.push_back(quad);

    } // for every vertex we have one quadrangle
  }

  require(triangles.size() * 3 == quadrangles.size(),
          "The number of quadrangles (" + d2s(quadrangles.size()) +
          " is not equal to number of triangles (" + d2s(triangles.size()) +
          " multiplying by 3 (" + d2s(3 * triangles.size()) + ")");

  // after that we check boundary element
  // because after adding new vertices they need to be redefined

  for (unsigned int line = 0; line < lines.size(); ++line)
  {
    // we need to find an edge that coincides with this line
    for (unsigned int edge = 0; edge < edges.size(); ++edge)
    {
      if (lines[line] == edges[edge])
      {
        // we change existing line and add new line in the end of list
        const unsigned int end_ver = lines[line].get_end();
        lines[line].set_end(n_old_vertices + edge); // change existing line
        // add new line
        lines.push_back(Line(n_old_vertices + edge,
                             end_ver,
                             lines[line].get_material_id()));
      }
    }
  }

}



void Mesh::edge_numeration()
{
  // matrix of incidence between vertices of the mesh
  const IncidenceMatrix incidence_matrix(vertices.size(), triangles);

  // the number of edges in such mesh - the number of non zero elements in incidence matrix
  const unsigned int n_edges = incidence_matrix.get_n_nonzero();

  // allocate memory for edges
  edges.resize(n_edges);

  // look through all cells of the mesh
  for (unsigned int cell = 0; cell < triangles.size(); ++cell)
  {
    unsigned int lne = 0; // local number of the edge (0 <= lne < cell::n_edges)
    for (unsigned int i = 0; i < Triangle::n_vertices; ++i)
    {
      const unsigned int ii = triangles[cell].get_vertex(i);
      for (unsigned int j = 0; j < Triangle::n_vertices; ++j)
      {
        const unsigned int jj = triangles[cell].get_vertex(j);
        if (ii > jj) // ii must be bigger than jj
        {
          const unsigned int gne = incidence_matrix.find(ii, jj); // global number of edge
          // set the global number of edge to cell
          triangles[cell].set_edge(lne, gne);
          // initialize edge
          edges[gne] = Line(std::min(ii, jj),
                            std::max(ii, jj),
                            triangles[cell].get_material_id());
          ++lne;
        }
      }
    }
    expect(lne == Triangle::n_edges,
           "lne must be equal to " + d2s(Triangle::n_edges) +
           ", but it is " + d2s(lne));
  }

} // edge numeration




void Mesh::write(const std::string &file)
{
  std::ofstream out(file.c_str());
  require(out, "File " + file + " cannot be opened for writing!");

  out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n" << vertices.size() << "\n";
  for (unsigned int ver = 0; ver < vertices.size(); ++ver)
  {
    out << ver + 1 << " ";
    for (unsigned int coord = 0; coord < Point::n_coord; ++coord)
      out << vertices[ver].get_coord(coord) << " ";
    out << "\n";
  }

//  const unsigned int n_elements = elements.size();
//  out << "$EndNodes\n$Elements\n" << n_elements << "\n";
//  for (unsigned int el = 0; el < n_elements; ++el)
//  {
//    out << el + 1                          /* serial number of element */
//        << elements[el].el_type_gmsh       /* type of element suitable for Gmsh */
//        << " 2 "                           /* the number of tags */
//        << elements[el].get_material_id()  /* physical domain */
//        << elements[el].get_material_id(); /* elemetary domain - let it be the same */
//    for (unsigned int ver = 0; ver < elements[el].n_vertices; ++ver)
//      out << elements[el].get_vertex(ver) << " ";
//    out << "\n";
//  }


  const unsigned int n_elements = quadrangles.size() + lines.size();
  out << "$EndNodes\n$Elements\n" << n_elements << "\n";
  for (unsigned int el = 0; el < quadrangles.size(); ++el)
  {
    out << el + 1 << " 3 2 " << quadrangles[el].get_material_id() << " 0 ";
    for (unsigned int ver = 0; ver < Quadrangle::n_vertices; ++ver)
      out << quadrangles[el].get_vertex(ver) + 1 << " ";
    out << "\n";
  }
  for (unsigned int el = 0; el < lines.size(); ++el)
  {
    out << quadrangles.size() + el + 1
        << " 1 2 " << lines[el].get_material_id()
        << " " << lines[el].get_material_id() << " "
        << lines[el].get_beg() + 1 << " " << lines[el].get_end() + 1;
    out << "\n";
  }

  out << "$EndElements\n";

  out.close();
}



inline unsigned int Mesh::get_n_vertices() const
{
  return vertices.size();
}

inline unsigned int Mesh::get_n_lines() const
{
  return lines.size();
}

inline unsigned int Mesh::get_n_edges() const
{
  return edges.size();
}

inline unsigned int Mesh::get_n_triangles() const
{
  return triangles.size();
}

inline unsigned int Mesh::get_n_tetrahedra() const
{
  return tetrahedra.size();
}

void Mesh::info(std::ostream &out) const
{
  out << "\nvertices    : " << vertices.size()
      << "\nedges       : " << edges.size()
      << "\nlines       : " << lines.size()
      << "\ntriangles   : " << triangles.size()
      << "\ntetrahedra  : " << tetrahedra.size()
      << "\nquadrangles : " << quadrangles.size()
      << "\n\n";
}


TETHEX_NAMESPACE_CLOSE
