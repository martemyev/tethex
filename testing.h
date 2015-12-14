#ifndef TETHEX_TESTING_H
#define TETHEX_TESTING_H

#include "config.h"

#if defined(TESTING)
#include "gtest/gtest.h"
#include "tethex.h"
#include <fstream>

using namespace tethex;


//---------------------------------------------------------
//
// Common auxiliary functions
//
//---------------------------------------------------------

void check_vertices(const Mesh &mesh,
                    const std::vector<Point> &vertices)
{
  EXPECT_EQ(mesh.get_n_vertices(), vertices.size());
  for (unsigned int i = 0; i < vertices.size(); ++i)
    for (unsigned int j = 0; j < Point::n_coord; ++j)
      EXPECT_DOUBLE_EQ(mesh.get_vertex(i).get_coord(j), vertices[i].get_coord(j));
}



void check_points(const Mesh &mesh,
                  const std::vector<PhysPoint> &points)
{
  EXPECT_EQ(mesh.get_n_points(), points.size());
  for (unsigned int i = 0; i < points.size(); ++i)
    for (unsigned int j = 0; j < PhysPoint::n_vertices; ++j)
      EXPECT_EQ(mesh.get_point(i).get_vertex(j), points[i].get_vertex(j));
}



void check_edges(const Mesh &mesh,
                 const std::vector<Line> &edges)
{
  EXPECT_EQ(mesh.get_n_edges(), edges.size());
  for (unsigned int i = 0; i < edges.size(); ++i)
    for (unsigned int j = 0; j < Line::n_vertices; ++j)
      EXPECT_EQ(mesh.get_edge(i).get_vertex(j), edges[i].get_vertex(j));
}



void check_lines(const Mesh &mesh,
                 const std::vector<Line> &lines)
{
  EXPECT_EQ(mesh.get_n_lines(), lines.size());
  for (unsigned int i = 0; i < lines.size(); ++i)
    for (unsigned int j = 0; j < Line::n_vertices; ++j)
      EXPECT_EQ(mesh.get_line(i).get_vertex(j), lines[i].get_vertex(j));
}



void check_faces(const Mesh &mesh,
                 const std::vector<Triangle> &faces)
{
  EXPECT_EQ(mesh.get_n_faces(), faces.size());
  for (unsigned int i = 0; i < faces.size(); ++i)
    for (unsigned int j = 0; j < Triangle::n_vertices; ++j)
    EXPECT_EQ(mesh.get_face(i).get_vertex(j), faces[i].get_vertex(j));
}



void check_triangle_vertices(const Mesh &mesh,
                             const std::vector<Triangle> &triangles)
{
  EXPECT_EQ(mesh.get_n_triangles(), triangles.size());
  for (unsigned int i = 0; i < triangles.size(); ++i)
    for (unsigned int j = 0; j < Triangle::n_vertices; ++j)
    EXPECT_EQ(mesh.get_triangle(i).get_vertex(j), triangles[i].get_vertex(j));
}



void check_triangle_edges(const Mesh &mesh,
                          const std::vector<Triangle> &triangles)
{
  EXPECT_EQ(mesh.get_n_triangles(), triangles.size());
  for (unsigned int i = 0; i < triangles.size(); ++i)
    for (unsigned int j = 0; j < Triangle::n_edges; ++j)
    EXPECT_EQ(mesh.get_triangle(i).get_edge(j), triangles[i].get_edge(j));
}



void check_quadrangle_vertices(const Mesh &mesh,
                               const std::vector<Quadrangle> &quadrangles)
{
  EXPECT_EQ(mesh.get_n_quadrangles(), quadrangles.size());
  for (unsigned int i = 0; i < quadrangles.size(); ++i)
    for (unsigned int j = 0; j < Quadrangle::n_vertices; ++j)
    EXPECT_EQ(mesh.get_quadrangle(i).get_vertex(j), quadrangles[i].get_vertex(j));
}



void check_tetrahedron_vertices(const Mesh &mesh,
                               const std::vector<Tetrahedron> &tetrahedra)
{
  EXPECT_EQ(mesh.get_n_tetrahedra(), tetrahedra.size());
  for (unsigned int i = 0; i < tetrahedra.size(); ++i)
    for (unsigned int j = 0; j < Tetrahedron::n_vertices; ++j)
    EXPECT_EQ(mesh.get_tetrahedron(i).get_vertex(j), tetrahedra[i].get_vertex(j));
}



void check_tetrahedron_edges(const Mesh &mesh,
                             const std::vector<Tetrahedron> &tetrahedra)
{
  EXPECT_EQ(mesh.get_n_tetrahedra(), tetrahedra.size());
  for (unsigned int i = 0; i < tetrahedra.size(); ++i)
    for (unsigned int j = 0; j < Tetrahedron::n_edges; ++j)
    EXPECT_EQ(mesh.get_tetrahedron(i).get_edge(j), tetrahedra[i].get_edge(j));
}



void check_tetrahedron_faces(const Mesh &mesh,
                             const std::vector<Tetrahedron> &tetrahedra)
{
  EXPECT_EQ(mesh.get_n_tetrahedra(), tetrahedra.size());
  for (unsigned int i = 0; i < tetrahedra.size(); ++i)
    for (unsigned int j = 0; j < Tetrahedron::n_faces; ++j)
    EXPECT_EQ(mesh.get_tetrahedron(i).get_face(j), tetrahedra[i].get_face(j));
}



void check_hexahedron_vertices(const Mesh &mesh,
                               const std::vector<Hexahedron> &hexahedra)
{
  EXPECT_EQ(mesh.get_n_hexahedra(), hexahedra.size());
  for (unsigned int i = 0; i < hexahedra.size(); ++i)
    for (unsigned int j = 0; j < Hexahedron::n_vertices; ++j)
    EXPECT_EQ(mesh.get_hexahedron(i).get_vertex(j), hexahedra[i].get_vertex(j));
}





//---------------------------------------------------------
//
// Common procedure 2D
//
//---------------------------------------------------------

void test_2D(const std::string &input_file,
             const std::string &output_file,
             std::vector<Point> &vertices,
             const std::vector<Point> &new_vertices,
             const std::vector<PhysPoint> points,
             const std::vector<Line> &edges,
             const std::vector<Line> &lines,
             const std::vector<Line> &new_lines,
             const std::vector<Triangle> &triangles,
             const std::vector<Quadrangle> &quadrangles)
{
  // read the mesh
  Mesh mesh;
  try
  {
    mesh.read(input_file);
  }
  catch (std::runtime_error &e)
  {
    std::cout << "\n=== EXCEPTION === \n";
    std::cout << e.what() << std::endl << std::endl;;
  }

  // check that mesh is right
  if (!vertices.empty())
    check_vertices(mesh, vertices);

  if (!lines.empty())
    check_lines(mesh, lines);

  if (!triangles.empty())
    check_triangle_vertices(mesh, triangles);

  // convert triangles into quadrangles
  try
  {
    mesh.convert();
  }
  catch (std::runtime_error &e)
  {
    std::cout << "\n=== EXCEPTION === \n";
    std::cout << e.what() << std::endl << std::endl;;
  }

  // check new vertices
  if (!new_vertices.empty())
  {
    vertices.insert(vertices.end(), new_vertices.begin(), new_vertices.end());
    check_vertices(mesh, vertices);
  }

  // during conversion we generate edges
  if (!edges.empty())
    check_edges(mesh, edges);

  // if we don't delete simplices,
  // we can compare edges of triangles
#if !defined(DELETE_SIMPLICES)
  if (!triangles.empty())
    check_triangle_edges(mesh, triangles);
#endif

  // now we compare quadrangles...
  if (!quadrangles.empty())
    check_quadrangle_vertices(mesh, quadrangles);

  // ... boundary lines ...
  if (!new_lines.empty())
    check_lines(mesh, new_lines);

  // ... and points
  if (!points.empty())
    check_points(mesh, points);

  // write the results for double checking
  if (!output_file.empty())
  {
    try
    {
      mesh.write(output_file);
    }
    catch (std::runtime_error &e)
    {
      std::cout << "\n=== EXCEPTION === \n";
      std::cout << e.what() << std::endl << std::endl;;
    }
  }
}



//---------------------------------------------------------
//
// Common procedure 3D
//
//---------------------------------------------------------

void test_3D(const std::string &input_file,
             const std::string &output_file,
             std::vector<Point> &vertices,
             const std::vector<Point> &new_vertices,
             const std::vector<PhysPoint> points,
             const std::vector<Line> &edges,
             const std::vector<Line> &lines,
             const std::vector<Line> &new_lines,
             const std::vector<Triangle> &faces,
             const std::vector<Triangle> &triangles,
             const std::vector<Quadrangle> &quadrangles,
             const std::vector<Tetrahedron> &tetrahedra,
             const std::vector<Hexahedron> &hexahedra)
{
  // read the mesh
  Mesh mesh;
  try
  {
    mesh.read(input_file);
  }
  catch (std::runtime_error &e)
  {
    std::cout << "\n=== EXCEPTION === \n";
    std::cout << e.what() << std::endl << std::endl;;
  }

  // check that mesh is right
  if (!vertices.empty())
    check_vertices(mesh, vertices);

  if (!lines.empty())
    check_lines(mesh, lines);

  if (!triangles.empty())
    check_triangle_vertices(mesh, triangles);

  if (!tetrahedra.empty())
    check_tetrahedron_vertices(mesh, tetrahedra);

  // convert tetrahedra into hexahedra
  try
  {
    mesh.convert();
  }
  catch (std::runtime_error &e)
  {
    std::cout << "\n=== EXCEPTION === \n";
    std::cout << e.what() << std::endl << std::endl;;
  }

  // check new vertices
  if (!new_vertices.empty())
  {
    vertices.insert(vertices.end(), new_vertices.begin(), new_vertices.end());
    check_vertices(mesh, vertices);
  }

  // during conversion we generate edges
  if (!edges.empty())
    check_edges(mesh, edges);

  // and faces
  if (!faces.empty())
    check_faces(mesh, faces);

  // if we don't delete simplices,
  // we can compare edges and faces of tetrahedra and
  // edges of boundary triangles
#if !defined(DELETE_SIMPLICES)
  if (!tetrahedra.empty())
  {
    check_tetrahedron_edges(mesh, tetrahedra);
    check_tetrahedron_faces(mesh, tetrahedra);
  }
#endif

  // now we compare hexahedra...
  if (!hexahedra.empty())
    check_hexahedron_vertices(mesh, hexahedra);

  // ... boundary quadrangles...
  if (!quadrangles.empty())
    check_quadrangle_vertices(mesh, quadrangles);

  // ... boundary lines ...
  if (!new_lines.empty())
    check_lines(mesh, new_lines);

  // ... and points
  if (!points.empty())
    check_points(mesh, points);

  // write the results for double checking
  if (!output_file.empty())
  {
    try
    {
      mesh.write(output_file);
    }
    catch (std::runtime_error &e)
    {
      std::cout << "\n=== EXCEPTION === \n";
      std::cout << e.what() << std::endl << std::endl;;
    }
  }
}





//---------------------------------------------------------
//
// Testing meshes declarations
//
//---------------------------------------------------------

void mesh_2d_1(std::string &input_file, std::string &output_file);
void mesh_2d_2(std::string &input_file, std::string &output_file);
void mesh_3d_1(std::string &input_file, std::string &output_file);
void mesh_3d_2(std::string &input_file, std::string &output_file);




//---------------------------------------------------------
//
// Tests 2D
//
//---------------------------------------------------------

TEST(Mesh_2D, OneElement)
{
  std::string input_file, output_file;
  mesh_2d_1(input_file, output_file);

  // vertices and mesh elements from the file
  std::vector<Point> vertices;
  vertices.push_back(Point(0, 0, 0));
  vertices.push_back(Point(1, 0, 0));
  vertices.push_back(Point(0, 1, 0));

  std::vector<Line> edges;
  edges.push_back(Line(0, 1, 11));
  edges.push_back(Line(0, 2, 11));
  edges.push_back(Line(1, 2, 11));

  std::vector<Line> lines;
  std::vector<Line> new_lines;

  std::vector<Triangle> triangles;
  triangles.push_back(Triangle(0, 1, 2, 11));
  triangles[0].set_edge(0, 0);
  triangles[0].set_edge(1, 1);
  triangles[0].set_edge(2, 2);

  std::vector<Quadrangle> quadrangles;
  quadrangles.push_back(Quadrangle(0, 3, 6, 4, 11));
  quadrangles.push_back(Quadrangle(1, 5, 6, 3, 11));
  quadrangles.push_back(Quadrangle(2, 4, 6, 5, 11));

  // new vertices created after conversion from triangles to quadrangles
  std::vector<Point> new_vertices;
  new_vertices.push_back(Point(1./2., 0, 0));
  new_vertices.push_back(Point(0, 1./2., 0));
  new_vertices.push_back(Point(1./2., 1./2., 0));
  new_vertices.push_back(Point(1./3., 1./3., 0));

  // physical points
  std::vector<PhysPoint> points;
  points.push_back(PhysPoint(0, 41));
  points.push_back(PhysPoint(2, 42));

  // testing
  test_2D(input_file,
          output_file,
          vertices,
          new_vertices,
          points,
          edges,
          lines,
          new_lines,
          triangles,
          quadrangles);

  // once again for converted mesh
  new_vertices.clear(); // vertices already include new_vertices
  edges.clear(); // we don't generate edges for quadrangles
  lines.clear(); // there are no original lines
#if defined(DELETE_SIMPLICES)
  triangles.clear(); // there are no triangles anymore
#endif
  test_2D(output_file,
          "",
          vertices,
          new_vertices,
          points,
          edges,
          lines,
          new_lines,
          triangles,
          quadrangles);
}




TEST(Mesh_2D, TwoElements)
{
  std::string input_file, output_file;
  mesh_2d_2(input_file, output_file);

  // vertices and mesh elements from the file
  std::vector<Point> vertices;
  vertices.push_back(Point(0, 0, 0));
  vertices.push_back(Point(1, 0, 0));
  vertices.push_back(Point(0, 1, 0));
  vertices.push_back(Point(1, 1, 0));

  std::vector<Line> edges;
  edges.push_back(Line(0, 1, 11));
  edges.push_back(Line(0, 2, 11));
  edges.push_back(Line(1, 2, 11));
  edges.push_back(Line(1, 3, 11));
  edges.push_back(Line(2, 3, 11));

  std::vector<Line> lines;
  lines.push_back(Line(0, 1, 21));
  lines.push_back(Line(0, 2, 22));
  lines.push_back(Line(1, 3, 23));
  lines.push_back(Line(2, 3, 24));

  std::vector<Line> new_lines;
  new_lines.push_back(Line(0, 4, 21));
  new_lines.push_back(Line(0, 5, 22));
  new_lines.push_back(Line(1, 7, 23));
  new_lines.push_back(Line(2, 8, 24));
  new_lines.push_back(Line(4, 1, 21));
  new_lines.push_back(Line(5, 2, 22));
  new_lines.push_back(Line(7, 3, 23));
  new_lines.push_back(Line(8, 3, 24));

  std::vector<Triangle> triangles;
  triangles.push_back(Triangle(0, 1, 2, 11));
  triangles[0].set_edge(0, 0);
  triangles[0].set_edge(1, 1);
  triangles[0].set_edge(2, 2);

  triangles.push_back(Triangle(1, 3, 2, 12));
  triangles[1].set_edge(0, 3);
  triangles[1].set_edge(1, 4);
  triangles[1].set_edge(2, 2);

  std::vector<Quadrangle> quadrangles;
  quadrangles.push_back(Quadrangle(0, 4, 9, 5, 11));
  quadrangles.push_back(Quadrangle(1, 6, 9, 4, 11));
  quadrangles.push_back(Quadrangle(2, 5, 9, 6, 11));
  quadrangles.push_back(Quadrangle(1, 7, 10, 6, 12));
  quadrangles.push_back(Quadrangle(3, 8, 10, 7, 12));
  quadrangles.push_back(Quadrangle(2, 6, 10, 8, 12));

  // new vertices created after conversion from triangles to quadrangles
  std::vector<Point> new_vertices;
  new_vertices.push_back(Point(1./2., 0, 0));
  new_vertices.push_back(Point(0, 1./2., 0));
  new_vertices.push_back(Point(1./2., 1./2., 0));
  new_vertices.push_back(Point(1, 1./2., 0));
  new_vertices.push_back(Point(1./2., 1, 0));
  new_vertices.push_back(Point(1./3., 1./3., 0));
  new_vertices.push_back(Point(2./3., 2./3., 0));

  // physical points
  std::vector<PhysPoint> points;

  // testing
  test_2D(input_file,
          output_file,
          vertices,
          new_vertices,
          points,
          edges,
          lines,
          new_lines,
          triangles,
          quadrangles);

  // once again for converted mesh
  new_vertices.clear(); // vertices already include new_vertices
  edges.clear(); // we don't generate edges for quadrangles
  lines.clear(); // there are no original lines
#if defined(DELETE_SIMPLICES)
  triangles.clear(); // there are no triangles anymore
#endif
  test_2D(output_file,
          "",
          vertices,
          new_vertices,
          points,
          edges,
          lines,
          new_lines,
          triangles,
          quadrangles);
}





//---------------------------------------------------------
//
// Tests 3D
//
//---------------------------------------------------------


TEST(Mesh_3D, OneElement)
{
  std::string input_file, output_file;
  mesh_3d_1(input_file, output_file);

  // vertices and mesh elements from the file
  std::vector<Point> vertices;
  vertices.push_back(Point(0, 0, 0));
  vertices.push_back(Point(1, 0, 0));
  vertices.push_back(Point(0, 1, 0));
  vertices.push_back(Point(0, 0, 1));

  // boundary elements are empty
  std::vector<Line> lines;
  std::vector<Line> new_lines;
  std::vector<Triangle> triangles;
  std::vector<Quadrangle> quadrangles;

  // mesh elements
  std::vector<Line> edges;
  edges.push_back(Line(0, 1, 11));
  edges.push_back(Line(0, 2, 11));
  edges.push_back(Line(1, 2, 11));
  edges.push_back(Line(0, 3, 11));
  edges.push_back(Line(1, 3, 11));
  edges.push_back(Line(2, 3, 11));

  std::vector<Triangle> faces;
  faces.push_back(Triangle(0, 1, 2, 11));
  faces.push_back(Triangle(0, 1, 3, 11));
  faces.push_back(Triangle(0, 2, 3, 11));
  faces.push_back(Triangle(1, 2, 3, 11));

  std::vector<Tetrahedron> tetrahedra;
  tetrahedra.push_back(Tetrahedron(0, 1, 2, 3, 11));
  tetrahedra[0].set_edge(0, 0);
  tetrahedra[0].set_edge(1, 1);
  tetrahedra[0].set_edge(2, 2);
  tetrahedra[0].set_edge(3, 3);
  tetrahedra[0].set_edge(4, 4);
  tetrahedra[0].set_edge(5, 5);
  tetrahedra[0].set_face(0, 0);
  tetrahedra[0].set_face(1, 1);
  tetrahedra[0].set_face(2, 2);
  tetrahedra[0].set_face(3, 3);

  // new vertices created after conversion from triangles to quadrangles
  std::vector<Point> new_vertices;
  new_vertices.push_back(Point(1./2., 0, 0));
  new_vertices.push_back(Point(0, 1./2., 0));
  new_vertices.push_back(Point(1./2., 1./2., 0));
  new_vertices.push_back(Point(0, 0, 1./2.));
  new_vertices.push_back(Point(1./2., 0, 1./2.));
  new_vertices.push_back(Point(0, 1./2., 1./2.));
  new_vertices.push_back(Point(1./3., 1./3., 0));
  new_vertices.push_back(Point(1./3., 0, 1./3.));
  new_vertices.push_back(Point(0, 1./3., 1./3.));
  new_vertices.push_back(Point(1./3., 1./3., 1./3.));
  new_vertices.push_back(Point(1./4., 1./4., 1./4.));

  std::vector<Hexahedron> hexahedra;
  hexahedra.push_back(Hexahedron(7, 11, 14, 12, 0, 4, 10, 5, 11));
  hexahedra.push_back(Hexahedron(1, 4, 10, 6, 8, 11, 14, 13, 11));
  hexahedra.push_back(Hexahedron(9, 12, 14, 13, 2, 5, 10, 6, 11));
  hexahedra.push_back(Hexahedron(3, 7, 11, 8, 9, 12, 14, 13, 11));

  // physical points
  std::vector<PhysPoint> points;
  points.push_back(PhysPoint(0, 41));
  points.push_back(PhysPoint(2, 42));
  points.push_back(PhysPoint(3, 43));

  // testing
  test_3D(input_file,
          output_file,
          vertices,
          new_vertices,
          points,
          edges,
          lines,
          new_lines,
          faces,
          triangles,
          quadrangles,
          tetrahedra,
          hexahedra);

  // once again for converted mesh
  new_vertices.clear(); // vertices already include new_vertices
  edges.clear(); // we don't generate edges for tetrahedra
  lines.clear(); // there are no original lines
  faces.clear(); // we don't generate faces for tetrahedra
#if defined(DELETE_SIMPLICES)
  triangles.clear(); // there are no triangles anymore
  tetrahedra.clear(); // and tetrahedra as well
#endif
  test_3D(output_file,
          "",
          vertices,
          new_vertices,
          points,
          edges,
          lines,
          new_lines,
          faces,
          triangles,
          quadrangles,
          tetrahedra,
          hexahedra);
}



TEST(Mesh_3D, TwoElements)
{
  std::string input_file, output_file;
  mesh_3d_2(input_file, output_file);

  // vertices and mesh elements from the file
  std::vector<Point> vertices;
  vertices.push_back(Point(0, 0, 0));
  vertices.push_back(Point(1, 0, 0));
  vertices.push_back(Point(0, 1, 0));
  vertices.push_back(Point(0, 0, 1));
  vertices.push_back(Point(0, 0, -1));

  std::vector<Line> edges;
  edges.push_back(Line(0, 1, 11));
  edges.push_back(Line(0, 2, 11));
  edges.push_back(Line(1, 2, 11));
  edges.push_back(Line(0, 3, 11));
  edges.push_back(Line(1, 3, 11));
  edges.push_back(Line(2, 3, 11));
  edges.push_back(Line(0, 4, 12));
  edges.push_back(Line(1, 4, 12));
  edges.push_back(Line(2, 4, 12));

  std::vector<Line> lines;
  lines.push_back(Line(0, 3, 31));
  lines.push_back(Line(4, 0, 32));

  std::vector<Line> new_lines;
  new_lines.push_back(Line(0, 8, 31));
  new_lines.push_back(Line(4, 11, 32));
  new_lines.push_back(Line(8, 3, 31));
  new_lines.push_back(Line(11, 0, 32));

  std::vector<Triangle> triangles;
  triangles.push_back(Triangle(0, 1, 3, 21));
  triangles.push_back(Triangle(4, 2, 0, 22));

  std::vector<Quadrangle> quadrangles;
  quadrangles.push_back(Quadrangle(0, 5, 15, 8, 21));
  quadrangles.push_back(Quadrangle(1, 5, 15, 9, 21));
  quadrangles.push_back(Quadrangle(3, 8, 15, 9, 21));
  quadrangles.push_back(Quadrangle(4, 13, 19, 11, 22));
  quadrangles.push_back(Quadrangle(2, 13, 19, 6, 22));
  quadrangles.push_back(Quadrangle(0, 11, 19, 6, 22));

  // new vertices created after conversion from triangles to quadrangles
  std::vector<Point> new_vertices;
  new_vertices.push_back(Point(1./2., 0, 0));
  new_vertices.push_back(Point(0, 1./2., 0));
  new_vertices.push_back(Point(1./2., 1./2., 0));
  new_vertices.push_back(Point(0, 0, 1./2.));
  new_vertices.push_back(Point(1./2., 0, 1./2.));
  new_vertices.push_back(Point(0, 1./2., 1./2.));
  new_vertices.push_back(Point(0, 0, -1./2.));
  new_vertices.push_back(Point(1./2., 0, -1./2.));
  new_vertices.push_back(Point(0, 1./2., -1./2.));
  new_vertices.push_back(Point(1./3., 1./3., 0));
  new_vertices.push_back(Point(1./3., 0, 1./3.));
  new_vertices.push_back(Point(0, 1./3., 1./3.));
  new_vertices.push_back(Point(1./3., 1./3., 1./3.));
  new_vertices.push_back(Point(1./3., 0, -1./3.));
  new_vertices.push_back(Point(0, 1./3., -1./3.));
  new_vertices.push_back(Point(1./3., 1./3., -1./3.));
  new_vertices.push_back(Point(1./4., 1./4., 1./4.));
  new_vertices.push_back(Point(1./4., 1./4., -1./4.));

  std::vector<Triangle> faces;
  faces.push_back(Triangle(0, 1, 2, 11));
  faces.push_back(Triangle(0, 1, 3, 11));
  faces.push_back(Triangle(0, 2, 3, 11));
  faces.push_back(Triangle(1, 2, 3, 11));
  faces.push_back(Triangle(0, 1, 4, 12));
  faces.push_back(Triangle(0, 2, 4, 12));
  faces.push_back(Triangle(1, 2, 4, 12));

  std::vector<Tetrahedron> tetrahedra;
  tetrahedra.push_back(Tetrahedron(0, 1, 2, 3, 11));
  tetrahedra[0].set_edge(0, 0);
  tetrahedra[0].set_edge(1, 1);
  tetrahedra[0].set_edge(2, 2);
  tetrahedra[0].set_edge(3, 3);
  tetrahedra[0].set_edge(4, 4);
  tetrahedra[0].set_edge(5, 5);
  tetrahedra[0].set_face(0, 0);
  tetrahedra[0].set_face(1, 1);
  tetrahedra[0].set_face(2, 2);
  tetrahedra[0].set_face(3, 3);

  tetrahedra.push_back(Tetrahedron(0, 1, 2, 4, 12));
  tetrahedra[1].set_edge(0, 0);
  tetrahedra[1].set_edge(1, 1);
  tetrahedra[1].set_edge(2, 2);
  tetrahedra[1].set_edge(3, 6);
  tetrahedra[1].set_edge(4, 7);
  tetrahedra[1].set_edge(5, 8);
  tetrahedra[1].set_face(0, 0);
  tetrahedra[1].set_face(1, 4);
  tetrahedra[1].set_face(2, 5);
  tetrahedra[1].set_face(3, 6);

  std::vector<Hexahedron> hexahedra;

  // physical points
  std::vector<PhysPoint> points;

  // testing
  test_3D(input_file,
          output_file,
          vertices,
          new_vertices,
          points,
          edges,
          lines,
          new_lines,
          faces,
          triangles,
          quadrangles,
          tetrahedra,
          hexahedra);

  // once again for converted mesh
  new_vertices.clear(); // vertices already include new_vertices
  edges.clear(); // we don't generate edges for tetrahedra
  lines.clear(); // there are no original lines
  faces.clear(); // we don't generate faces for tetrahedra
#if defined(DELETE_SIMPLICES)
  triangles.clear(); // there are no triangles anymore
  tetrahedra.clear(); // and tetrahedra as well
#endif
  test_3D(output_file,
          "",
          vertices,
          new_vertices,
          points,
          edges,
          lines,
          new_lines,
          faces,
          triangles,
          quadrangles,
          tetrahedra,
          hexahedra);
}



//---------------------------------------------------------
//
// Testing meshes definitions
//
//---------------------------------------------------------

void mesh_2d_1(std::string &input_file, std::string &output_file)
{
  input_file = "mesh_2d_1.msh";
  output_file = "mesh_2d_1_hex.msh";
  std::ofstream out(input_file.c_str());
  require(out, "File " + input_file + " cannot be opened!");
  out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
  out << "3\n";
  out << "1 0 0 0\n";
  out << "2 1 0 0\n";
  out << "3 0 1 0\n";
  out << "$EndNodes\n$Elements\n";
  out << "3\n";
  out << "1 2 2 11 0 1 2 3\n";
  out << "2 15 2 41 0 1\n";
  out << "3 15 2 42 0 3\n";
  out << "$EndElements\n";
  out.close();
}




void mesh_2d_2(std::string &input_file, std::string &output_file)
{
  input_file = "mesh_2d_2.msh";
  output_file = "mesh_2d_2_hex.msh";
  std::ofstream out(input_file.c_str());
  require(out, "File " + input_file + " cannot be opened!");
  out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
  out << "4\n";
  out << "1 0 0 0\n";
  out << "2 1 0 0\n";
  out << "3 0 1 0\n";
  out << "4 1 1 0\n";
  out << "$EndNodes\n$Elements\n";
  out << "6\n";
  out << "1 2 2 11 0 1 2 3\n";
  out << "2 2 2 12 0 2 4 3\n";
  out << "3 1 2 21 0 1 2\n";
  out << "4 1 2 22 0 1 3\n";
  out << "5 1 2 23 0 2 4\n";
  out << "6 1 2 24 0 3 4\n";
  //out << "3 15 2 21 0 1\n"; // physical point
  out << "$EndElements\n";
  out.close();
}




void mesh_3d_1(std::string &input_file, std::string &output_file)
{
  input_file = "mesh_3d_1.msh";
  output_file = "mesh_3d_1_hex.msh";
  std::ofstream out(input_file.c_str());
  require(out, "File " + input_file + " cannot be opened!");
  out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
  out << "4\n";
  out << "1 0 0 0\n";
  out << "2 1 0 0\n";
  out << "3 0 1 0\n";
  out << "4 0 0 1\n";
  out << "$EndNodes\n$Elements\n";
  out << "4\n";
  out << "1 4 2 11 0 1 2 3 4\n";
  out << "2 15 2 41 0 1\n";
  out << "3 15 2 42 0 3\n";
  out << "4 15 2 43 0 4\n";
  out << "$EndElements\n";
  out.close();
}




void mesh_3d_2(std::string &input_file, std::string &output_file)
{
  input_file = "mesh_3d_2.msh";
  output_file = "mesh_3d_2_hex.msh";
  std::ofstream out(input_file.c_str());
  require(out, "File " + input_file + " cannot be opened!");
  out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
  out << "5\n";
  out << "1 0 0 0\n";
  out << "2 1 0 0\n";
  out << "3 0 1 0\n";
  out << "4 0 0 1\n";
  out << "5 0 0 -1\n";
  out << "$EndNodes\n$Elements\n";
  out << "6\n";
  out << "1 4 2 11 0 1 2 3 4\n";
  out << "2 4 2 12 0 1 2 3 5\n";
  out << "3 2 2 21 0 1 2 4\n";
  out << "4 2 2 22 0 5 3 1\n";
  out << "5 1 2 31 0 1 4\n";
  out << "6 1 2 32 0 5 1\n";
  out << "$EndElements\n";
  out.close();
}

#endif // TESTING

#endif // TETHEX_TESTING_H
