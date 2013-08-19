cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type of the project (Debug, Release)")


# --- some definitions ---
if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-DDEBUG)
endif(CMAKE_BUILD_TYPE MATCHES Debug)


# --- project ---
project(tethex)


# --- set version ---
set(TETHEX_VERSION_MAJOR 0)
set(TETHEX_VERSION_MINOR 0)
set(TETHEX_VERSION_PATCH 4)


# --- useful variables and options ---
option(TESTING "Do you want to launch testing procedure before using the program?" OFF)
option(DELETE_SIMPLICES "Delete simplices after converting them into bricks (to economize some memory)" ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

set(HAVE_64BIT_SIZE_T OFF CACHE INTERNAL "")
include(CheckTypeSize)
check_type_size("void*" SIZEOF_VOID_P)
if(SIZEOF_VOID_P EQUAL 8)
  set(HAVE_64BIT_SIZE_T ON CACHE INTERNAL "")
endif(SIZEOF_VOID_P EQUAL 8)


# --- write the variables initialized here to config.h file ---
message("Options:")
message("  CMAKE_BUILD_TYPE  = ${CMAKE_BUILD_TYPE}")
message("  TESTING           = ${TESTING}")
message("  DELETE_SIMPLICES  = ${DELETE_SIMPLICES}")
message("  HAVE_64BIT_SIZE_T = ${HAVE_64BIT_SIZE_T}")
configure_file(
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h")
include_directories(${PROJECT_BINARY_DIR}) # because config.h file is in binary directory we need to include it


# --- sources and headers ---
set(SRC_LIST main.cpp tethex.cpp)
set(HDR_LIST tethex.h)


# --- Google Test ---
if(TESTING)
  message("ATTENTION! This is a version with testing features!")
  set(GTEST_DIR /home/artemiev/projects/gtest-1.6.0 CACHE PATH "Path to Google Test top level directory")
  find_library(GTEST_LIB gtest ${GTEST_DIR}/lib/.libs)
  if(GTEST_LIB)
    message("GTest library is ${GTEST_LIB}")
  else(GTEST_LIB)
    message(FATAL_ERROR "GTest was not found!")
  endif(GTEST_LIB)
  include_directories(${GTEST_DIR}/include)
endif(TESTING)


# --- make ---
add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})
if(TESTING)
  target_link_libraries(${PROJECT_NAME} ${GTEST_LIB} pthread)
endif(TESTING)


# --- packaging of the library ---
set(CPACK_PACKAGE_VERSION_MAJOR "${TETHEX_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${TETHEX_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${TETHEX_VERSION_PATCH}")
# -- make package_source --
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "./*.user;build;doc;.git;.gitignore;dont_remove")
include(CPack)


# --- autodocumentation using doxygen ---
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(
    "${PROJECT_SOURCE_DIR}/doxyfile.in"
    "${CMAKE_CURRENT_BINARY_DIR}/doxyfile")
  add_custom_target(doc COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxyfile)
else(DOXYGEN_FOUND)
  message("WARNING: Doxygen is not found - Autodocumentation will not be created")
endif(DOXYGEN_FOUND)


# --- ctest and cdash ---
include(CTest)

