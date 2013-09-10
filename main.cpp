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
#include <ctime>
#if defined(TESTING)
  #include "testing.h"
#endif
#if defined(BOOST)
  #include <boost/program_options.hpp>
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

  run(in_file, out_file, verbosity);

  return 0;
}



int read_parameters(int argc, char **argv,
                    std::string &in_file,
                    std::string &out_file,
                    int &verbosity)
{
#if defined(BOOST)
  namespace po = boost::program_options;

  // declare the supported options
  po::options_description desc("\nProgram usage: ./tethex [arguments]\nHere is a list of allowed arguments");
  desc.add_options()
      ("input,i", po::value<std::string>(),
       ("name of input file (mandatory)"))
      ("output,o", po::value<std::string>(),
       ("name of output file (optional; if omitted, input file name is used with '_hex' before extension)"))
      ("verbosity", po::value<int>(),
       ("set verbosity level: 0 - no messages, 1 - short mesh statistics (by default), 2 - detailed mesh statistics (experimental)"))
      ("help,h", "produce help message")
      ("version,v", "print version")
      ("info", "print detailed information about the program")
      ("test", "launch testing before the conversion")
  ;

  po::variables_map vm;
  po::store(
        po::command_line_parser(argc, argv)
          .options(desc)
          .style(
              (po::command_line_style::unix_style |
               po::command_line_style::allow_long_disguise)
             ^(po::command_line_style::allow_guessing))
          .run(),
        vm);
  po::notify(vm);

  if (argc == 1 || vm.count("help")) // there is no arguments at all or there is a request for help
  {
    std::cout << desc << std::endl;
    return 1;
  }

  if (vm.count("version")) // request for version
  {
    std::cout << TETHEX_VERSION << std::endl;
    return 1;
  }

  if (vm.count("info")) // request for information about the program
  {
    std::string info = "tethex " + std::string(TETHEX_VERSION) + "\n";
    info += "build type: ";
#if defined(DEBUG)
    info += "Debug\n";
#else
    info += "Release\n";
#endif
    info += "options: ";
#if defined(DEBUG)
    info += "DEBUG ";
#endif
#if defined(TESTING)
    info += "TESTING ";
#endif
#if defined(DELETE_SIMPLICES)
    info += "DELETE_SIMPLICES ";
#endif
#if defined(HAVE_64BIT_SIZE_T)
    info += "HAVE_64BIT_SIZE_T ";
#endif
    info += "BOOST ";
    std::cout << info << std::endl;
    return 1;
  }

  if (vm.count("test")) // include testing in the conversion process
  {
#if defined(TESTING)
    std::cout << "\nWe are starting short testing procedures!\n\n";
    ::testing::InitGoogleTest(&argc, argv);
    int test_ret = RUN_ALL_TESTS();
    std::cout << "\nTesting procedures finished (" << test_ret << " is returned)\n\n";
#else
    std::cout << "\ntethex has been built without testing support.\n" << std::endl;
#endif
  }

  unsigned int extension_pos;
  if (vm.count("input"))
  {
    in_file = vm["input"].as<std::string>();
    extension_pos = in_file.find(".msh");
    require(extension_pos != in_file.npos,
            "Input file " + in_file +
            " is not native Gmsh's mesh file. Its extension is not .msh");
  }
  else // this parameter is mandatory
    require(false, "There is no input mesh file. Please use '-help' option to check how you can define an input mesh file");

  if (vm.count("output"))
    out_file = vm["output"].as<std::string>();
  else
    out_file = in_file.replace(extension_pos, in_file.size(), "_hex.msh");

  verbosity = 1; // by default
  if (vm.count("verbosity"))
  {
    verbosity = vm["verbosity"].as<int>();
    require(verbosity >= 0 && verbosity <= 2, "Incorrect verbosity value (" + d2s(verbosity) + ")");
  }

#else // BOOST

#if defined(TESTING)
  std::cout << "\nWe are starting short testing procedures!\n\n";
  ::testing::InitGoogleTest(&argc, argv);
  int test_ret = RUN_ALL_TESTS();
  std::cout << "\nTesting procedures finished (" << test_ret << " is returned)\n\n";
#endif

  if (argc < 2)
  {
    std::cout << "\nThere must be at least one argument - the name of mesh file." << std::endl;
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

#endif // BOOST

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


