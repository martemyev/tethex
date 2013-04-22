#ifndef TETHEX_H
#define TETHEX_H

#include <string>
#include <stdexcept>
#include <sstream>
#include <vector>

#define TETHEX_NAMESPACE_OPEN namespace tethex {
#define TETHEX_NAMESPACE_CLOSE }

TETHEX_NAMESPACE_OPEN

//-------------------------------------------------------
//
// d2s - convert data to string
//
//-------------------------------------------------------
                /**
                 * Convert float number to string
                 * @param x - the number in double format
                 * @param scientific - use scientific format (e.g. 1e-10), or not
                 * @param precision - if scientific format is used, we can change the precision
                 * @return string
                 */
std::string d2s(double x, bool scientific = false, int precision = 6);

                /**
                 * Convert integer data to string
                 * @param x - the integer number
                 * @return string
                 */
std::string d2s(int x);

                /**
                 * Convert unsigned integer data to string
                 * @param x - unsigned integer number
                 * @return string
                 */
std::string d2s(unsigned int x);




//-------------------------------------------------------
//
// expect and require
//
//-------------------------------------------------------
#if DEBUG
  #define expect(condition, message) \
    if (!(condition))                \
      requirement_fails(__FILE__,     \
                        __LINE__,     \
                        message)
#else
  // in release (or release-like) versions
  // nothing happens
  #define expect(condition, message) { }
#endif // DEBUG

#define require(condition, message) \
  if (!(condition))                 \
    requirement_fails(__FILE__,     \
                      __LINE__,     \
                      message)

                /**
                 * Throw an informative exception,
                 * if requirement or expectation fails
                 */
void requirement_fails(const char *file,
                       unsigned int line,
                       std::string message);




//-------------------------------------------------------
//
// Point
//
//-------------------------------------------------------
/**
 * Point in 3-dimensional space.
 * It's used for 2-dimensional triangles as well, and
 * in this case third coordinate is 0.
 */
class Point
{
public:
                /**
                 * The number of Cartesian coordinates, that describe the point.
                 * Here we always use 3 coordinates to describe a point.
                 */
  static const unsigned int n_coord = 3;

                /**
                 * Constructor
                 */
  Point();

                /**
                 * Constructor
                 */
  Point(const double coordinates[]);

                /**
                 * Copy constructor
                 */
  Point(const Point &p);

                /**
                 * Copy assignment operator
                 */
  Point& operator =(const Point &p);

                /**
                 * Get the coordinate of the point
                 * @param number - the serial number of coordinate [0, n_coord)
                 */
  double get_coord(unsigned int number) const;

                /**
                 * Set the value of specific coordinate
                 * @param number - the number of coordinate that we want to set
                 * @param value - new value of coordinate
                 */
  void set_coord(unsigned int number, double value);

private:
                /**
                 * Cartesian coordinates of the point
                 */
  double coord[n_coord];
};





//-------------------------------------------------------
//
// MeshElement
//
//-------------------------------------------------------
class MeshElement
{
public:
                /**
                 * Get the number of vertices
                 */
  unsigned int get_n_vertices() const;

                /**
                 * Get the number of edges
                 */
  unsigned int get_n_edges() const;

                /**
                 * Get type of the element that is used in Gmsh
                 */
  unsigned int get_gmsh_el_type() const;

                /**
                 * Get the material ID of the element
                 * @return The number that describes the physical domain to which the element belongs
                 */
  unsigned int get_material_id() const;

                /**
                 * Get the number of vertex describing the element
                 * @param number - local number of vertex [0, n_vertices)
                 * @return global number of vertex (among other mesh vertices)
                 */
  unsigned int get_vertex(unsigned int number) const;

                /**
                 * Get the number of edge describing the element
                 * @param number - local number of edge [0, n_edges)
                 * @return global number of edge (among other mesh edges)
                 */
  unsigned int get_edge(unsigned int number) const;

                /**
                 * Set the number of vertex
                 * @param local_number - the number of vertex inside the element [0, n_vertices)
                 * @param global_unmber - the number of vertex among other vertices of the mesh
                 */
  void set_vertex(unsigned int local_number, unsigned int global_number);

                /**
                 * Set the number of edge
                 * @param local_number - the number of edge inside the element [0, n_edges)
                 * @param global_unmber - the number of edge among other edges of the mesh
                 */
  void set_edge(unsigned int local_number, unsigned int global_number);

protected:
                /**
                 * The number of vertices describing the element.
                 * It must be defined in each derived class,
                 * because it's 0 by default.
                 */
  unsigned int n_vertices;

                /**
                 * Vertices (i.e. their global numbers) describing the element
                 */
  std::vector<unsigned int> vertices;

                /**
                 * The number of edges describing the element.
                 * It must be defined in each derived class,
                 * because it's 0 by default.
                 */
  unsigned int n_edges;

                /**
                 * Edges (i.e. their global numbers) describing the element
                 * It's not always used.
                 */
  std::vector<unsigned int> edges;

                /**
                 * ID of the physical domain where the element takes place.
                 * It's necessary to distinguish media with different physical properties.
                 */
  unsigned int material_id;

                /**
                 * Type of the element (its number actually) like in Gmsh.
                 * It must be defined in every derived class.
                 * It's 0 be default.
                 */
  unsigned int gmsh_el_type;

                /**
                 * Constructor is protected to prevent creating MeshElement objects directly
                 * @param n_ver - number of vertices
                 * @param n_edg - number of edges
                 * @param el_type - type of the element in Gmsh
                 */
  MeshElement(unsigned int n_ver = 0,
              unsigned int n_edg = 0,
              unsigned int el_type = 0);

                /**
                 * Copy constructor
                 */
  MeshElement(const MeshElement &elem);

                /**
                 * Copy assignment operator
                 */
  MeshElement& operator =(const MeshElement &elem);

                /**
                 * Destructor
                 */
  virtual ~MeshElement();
};




//-------------------------------------------------------
//
// Line
//
//-------------------------------------------------------
/**
 * Line keeps 2 numbers - the number of beginning vertex
 * and the number of the ending one, and a number of physical domain
 * where this line takes place (material identificator - in other words).
 * Therefore, Line is like an Edge but vertices are not ordered, and
 * plus material_id exists.
 */
class Line : public MeshElement
{
public:
                /**
                 * There are 2 vertices to describe a line
                 */
  static const unsigned int n_vertices = 2;

                /**
                 * Line is edge itself, so the number of edges is 1
                 */
  static const unsigned int n_edges = 1;

                /**
                 * In Gmsh line (physical line) is defined by number 1
                 */
  static const unsigned int gmsh_el_type = 1;

                /**
                 * Constructor
                 */
  Line();

                /**
                 * Constructor with parameters
                 * @param ver - the list of vertices
                 * @param mat_id - the material ID
                 */
  Line(const std::vector<unsigned int> &ver,
       const unsigned int mat_id = 0);

                /**
                 * Constructor with parameters
                 * @param v1 - one vertex
                 * @param v2 - another vertex
                 * @param mat_id - material ID
                 */
  Line(const unsigned int v1,
       const unsigned int v2,
       const unsigned int mat_id = 0);

                /**
                 * Comparing lines by their vertices
                 */
  bool operator ==(const Line &line) const;
};




//-------------------------------------------------------
//
// Triangle
//
//-------------------------------------------------------
/**
 * Triangle - 2-dimensional simplex
 */
class Triangle : public MeshElement
{
public:
                /**
                 * The number of vertices of triangle
                 */
  static const unsigned int n_vertices = 3;

                /**
                 * The number of edges of triangle
                 */
  static const unsigned int n_edges = 3;

                /**
                 * In Gmsh triangle is defined by number 2
                 */
  static const unsigned int gmsh_el_type = 2;

                /**
                 * Default constructor
                 */
  Triangle();

                /**
                 * Constructor with parameters
                 * @param ver - triangle vertices
                 * @param mat_id - material ID
                 */
  Triangle(const std::vector<unsigned int> &ver,
           const unsigned int mat_id = 0);
};




//-------------------------------------------------------
//
// Tetrahedron
//
//-------------------------------------------------------
/**
 * Tetrahedron - 3-dimensional simplex
 */
class Tetrahedron : public MeshElement
{
public:
                /**
                 * the number of vertices of tetrahedron
                 */
  static const unsigned int n_vertices = 4;

                /**
                 * the number of edges of tetrahedron
                 */
  static const unsigned int n_edges = 6;

                /**
                 * In Gmsh triangle is defined by number 2
                 */
  static const unsigned int gmsh_el_type = 4;

                /**
                 * Default constructor.
                 */
  Tetrahedron();

                /**
                 * Constructor with parameters
                 * @param ver - triangle vertices
                 * @param mat_id - material ID
                 */
  Tetrahedron(const std::vector<unsigned int> &ver,
              const unsigned int mat_id = 0);
};




//-------------------------------------------------------
//
// Quadrangle
//
//-------------------------------------------------------
/**
 * Quadrangle - 2-dimensional shape with 4 straight edges
 */
class Quadrangle : public MeshElement
{
public:
                /**
                 * the number of vertices of quadrangle
                 */
  static const unsigned int n_vertices = 4;

                /**
                 * the number of edges of quadrangle
                 */
  static const unsigned int n_edges = 4;

                /**
                 * In Gmsh triangle is defined by number 2
                 */
  static const unsigned int gmsh_el_type = 3;

                /**
                 * Default constructor.
                 */
  Quadrangle();

                /**
                 * Constructor with parameters
                 * @param ver - triangle vertices
                 * @param mat_id - material ID
                 */
  Quadrangle(const std::vector<unsigned int> &ver,
             const unsigned int mat_id = 0);
};







//-------------------------------------------------------
//
// IncidenceMatrix
//
//-------------------------------------------------------

//class MeshElement;

/**
 * Incidence matrix describes the relations (connections) between mesh nodes.
 * It's used for edge numerarion.
 * Because this matrix is symmetric we consider its lower triangle only.
 */
class IncidenceMatrix
{
public:
                /**
                 * Constructor
                 * @param n_vertices - the number of all mesh vertices
                 * @param cells - the list of all mesh cells
                 */
  IncidenceMatrix(const unsigned int n_vertices,
                  const std::vector<Triangle> &cells);

                /**
                 * Destructor
                 */
  ~IncidenceMatrix();

                /**
                 * Find a serial number of the non zero element in the matrix
                 * @param row_number - the number of row where we seek
                 * @param col_number - the number of column where we seek
                 * @return Serial number of the non zero element in the matrix
                 */
  unsigned int find(const unsigned int row_number,
                    const unsigned int col_number) const;

                /**
                 * Get the number of non zero elements in the matrix
                 */
  unsigned int get_n_nonzero() const;

private:
                /**
                 * The dimension of the matrix
                 */
  unsigned int dim;

                /**
                 * The number of nonzero elements of lower matrix triangle
                 */
  unsigned int n_non_zero;

                /**
                 * The number of nonzero elements in each row of lower matrix triangle
                 */
  unsigned int *row;

                /**
                 * The numbers of nonzero elements of lower matrix triangle
                 */
  unsigned int *col;
};




//-------------------------------------------------------
//
// Mesh
//
//-------------------------------------------------------
class Mesh
{
public:
                /**
                 * Constructor - nothing special
                 */
  Mesh();

                /**
                 * Destructor - to clean the memory
                 */
  ~Mesh();

                /**
                 * Read the mesh from file
                 * @param file - the name of the mesh file
                 */
  void read(const std::string &file);

                /**
                 * Conversion from simplices to bricks
                 */
  void convert();

                /**
                 * Write the resulting brick mesh to file
                 * @param file - the name of mesh file where we write the results of conversion
                 */
  void write(const std::string &file);

                /**
                 * Get the number of vertices
                 */
  unsigned int get_n_vertices() const;

                /**
                 * Get the number of lines
                 */
  unsigned int get_n_lines() const;

                /**
                 * Get the number of edges
                 */
  unsigned int get_n_edges() const;

                /**
                 * Get the number of triangles
                 */
  unsigned int get_n_triangles() const;

                /**
                 * Get the number of tetrahedra
                 */
  unsigned int get_n_tetrahedra() const;

  //unsigned int get_n_elements() const;

                /**
                 * Print information about mesh
                 */
  void info(std::ostream &out) const;

private:
                /**
                 * Mesh vertices (nodes)
                 */
  std::vector<Point> vertices;

                /**
                 * Mesh lines - mean physical lines
                 */
  std::vector<Line> lines;

                /**
                 * Mesh edges (oriented lines)
                 */
  //std::vector<Edge> edges;
  std::vector<Line> edges;

                /**
                 * Mesh triangles
                 */
  std::vector<Triangle> triangles;

                /**
                 * Mesh tetrahedra
                 */
  std::vector<Tetrahedron> tetrahedra;

                /**
                 * Mesh quadrangles
                 */
  std::vector<Quadrangle> quadrangles;

                /**
                 * Free the memory to read again, for example
                 */
  void clean();

                /**
                 * Numerate the edges
                 */
  void edge_numeration();

};


TETHEX_NAMESPACE_CLOSE

#endif // TETHEX_H
