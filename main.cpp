/*
 * tethex - tetrahedra to hexahedra conversion
 * Copyright (c) 2013 Mikhail Artemiev
 *
 * http://code.google.com/p/tethex
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
#include <iostream>
#include <cstdlib>
#if defined(TESTING)
  #include "testing.h"
#endif

int main(int argc, char **argv)
{
#if defined(TESTING)
  std::cout << "\nWe are starting short testing procedures!\n\n";
  ::testing::InitGoogleTest(&argc, argv);
  int test_ret = RUN_ALL_TESTS();
  std::cout << "\nTesting procedures finished (" << test_ret << " is returned)\n\n";
#endif

  using namespace tethex;

  if (argc < 2)
  {
    std::cout << "\nThere must be at least one argument - the name of mesh file." << std::endl;
    std::cout << "example of using: ./tethex input.msh [output.msh]" << std::endl;
    std::cout << "output.msh is the name of file containing converted elements." << std::endl;
    std::cout << "If output.msh is omitted, the resulting file will be named input_hex.msh.\n" <<std::endl;
    return 0;
  }

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

  return 0;
}

