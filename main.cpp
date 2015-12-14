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
#include <iostream>
#include <cstdlib>
#include <ctime>
#if defined(TESTING)
  #include "testing.h"
#endif

using namespace tethex;

int read_parameters(int argc, char **argv,
                    std::string &in_file,
                    std::string &out_file,
                    int &verbosity);

void run(const std::string &in_file,
         const std::string &out_file,
         int verbosity);



int main(int argc, char **argv)
{
  std::string in_file;  // input file
  std::string out_file; // output file
  int verbosity;        // verbosity level

  if(read_parameters(argc, argv, in_file, out_file, verbosity))
    return 1;

#if defined(TESTING)
  std::cout << "\nWe are starting short testing procedures!\n\n";
  ::testing::InitGoogleTest(&argc, argv);
  int test_ret = RUN_ALL_TESTS();
  std::cout << "\nTesting procedures finished (" << test_ret << " is returned)\n\n";
#endif

  run(in_file, out_file, verbosity);

  return 0;
}



int read_parameters(int argc, char **argv,
                    std::string &in_file,
                    std::string &out_file,
                    int &verbosity)
{
  if (argc < 2)
  {
    std::cout << "\nThere must be at least one argument - the name of a mesh file." << std::endl;
    std::cout << "example of using: ./tethex input.msh [output.msh]" << std::endl;
    std::cout << "output.msh is the name of file containing converted elements." << std::endl;
    std::cout << "If output.msh is omitted, the resulting file will be named input_hex.msh.\n" << std::endl;
    return 1;
  }

  in_file = std::string(argv[1]);
  const unsigned int pos = in_file.find(".msh");
  require(pos != in_file.npos,
          "Input file " + in_file +
          " is not native Gmsh's mesh file. Its extension is not .msh");

  if (argc > 2)
    out_file = std::string(argv[2]);
  else
    out_file = in_file.replace(pos, in_file.size(), "_hex.msh");

  verbosity = 1;

  return 0;
}



void run(const std::string &in_file,
         const std::string &out_file,
         int verbosity)
{
  // time measurement
  clock_t beg_time = clock();

  Mesh mesh;
  if (verbosity)
    std::cout << "Reading " << in_file << " file..." << std::endl;
  mesh.read(in_file);
  if (verbosity)
    std::cout << "Reading " << in_file << " file is done" << std::endl;

  if (verbosity) mesh.info();
  if (verbosity > 1) mesh.statistics(std::cout);

  if (verbosity)
    std::cout << "Converting simplices to bricks..." << std::endl;
  mesh.convert();
  if (verbosity)
    std::cout << "Converting simplices to bricks is done" << std::endl;

  if (verbosity) mesh.info();
  if (verbosity > 1) mesh.statistics(std::cout);

  if (verbosity)
    std::cout << "Writing " << out_file << " file..." << std::endl;
  mesh.write(out_file);
  if (verbosity)
    std::cout << "Writing " << out_file << " file is done" << std::endl;

  // time measurement
  clock_t total_time = clock() - beg_time;
  if (verbosity)
    std::cout << "\nTime of execution is " << (double)total_time / CLOCKS_PER_SEC << " sec " << std::endl;
}


