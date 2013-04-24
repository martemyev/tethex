#include "tethex.h"
#include "config.h"
#include <iostream>
#include <cstdlib>
#if TESTING
  #include "testing.h"
#endif

int main(int argc, char **argv)
{
#if TESTING
  std::cout << "\nWe are starting short testing procedures!\n\n";
  ::testing::InitGoogleTest(&argc, argv);
  int test_ret = RUN_ALL_TESTS();
  std::cout << "\nThe returning values of all tests running is " << test_ret << "\n\n";
#endif

  using namespace tethex;

  require(argc > 1,
          "There must be at least one argument with the name of mesh file");
  std::string file_name = std::string(argv[1]);
  const unsigned int pos = file_name.find(".msh");
  require(pos != std::string::npos,
          "Input file " + file_name +
          " is not native Gmsh's mesh file. Its extension is not .msh");

  Mesh mesh;
  std::cout << "Reading " << file_name << " file..." << std::endl;
  mesh.read(file_name);
  std::cout << "Reading " << file_name << " file is done" << std::endl;

  mesh.info(std::cout);
//  mesh.statistics(std::cout);

  std::cout << "Converting simplices to bricks..." << std::endl;
  mesh.convert();
  std::cout << "Converting simplices to bricks is done" << std::endl;

  mesh.info(std::cout);
//  mesh.statistics(std::cout);

  std::string res_name = file_name.replace(pos, file_name.size(), "_hex.msh");
  std::cout << "Writing " << res_name << " file..." << std::endl;
  mesh.write(res_name);
  std::cout << "Writing " << res_name << " file is done" << std::endl;

  std::string tmp = "optirun gmsh " + res_name;
  system(tmp.c_str());

  return 0;
}

