cmake_minimum_required(VERSION 2.8)


# --- some definitions ---
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type of the project (Debug, Release)")
set(DEBUG OFF CACHE BOOL "Build the library with debug features! Don't confuse with compiler debug options!")

if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(DEBUG ON CACHE BOOL "Build the library with debug features! Don't confuse with compiler debug options!" FORCE)
endif()


# --- project ---
project(tethex CXX)


# --- set version ---
set(TETHEX_VERSION_MAJOR 0)
set(TETHEX_VERSION_MINOR 0)
set(TETHEX_VERSION_PATCH 7)


# --- useful variables and options ---
option(TESTING "Do you want to launch testing procedure before using the program?" OFF)
option(DELETE_SIMPLICES "Delete simplices after converting them into bricks (to save some memory)" ON)
if(TESTING AND NOT DELETE_SIMPLICES)
  message("If TESTING=ON we need to remove simplices. Otherwise testing procedures will fail. Therefore DELETE_SIMPLICES=ON")
  set(DELETE_SIMPLICES ON CACHE BOOL "" FORCE)
endif()


# --- C++11 ---
if(CMAKE_COMPILER_IS_GNUCXX)
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
  if(GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
    message(STATUS "C++11 is activated")
    add_definitions("-std=gnu++11")
  elseif(GCC_VERSION VERSION_GREATER 4.3 OR GCC_VERSION VERSION_EQUAL 4.3)
    message(STATUS "C++0x is activated")
    add_definitions("-std=gnu++0x")
  else()
    message(FATAL_ERROR "C++11 or C++0x needed. Therefore a gcc compiler with a version higher than (or equal to) 4.3 is needed.")
  endif()
else()
  message(STATUS "C++0x is activated")
  add_definitions("-std=c++0x")
endif()


# --- write the variables initialized here to config.h file ---
message("Options:")
message("  CMAKE_BUILD_TYPE  = ${CMAKE_BUILD_TYPE}")
message("  TESTING           = ${TESTING}")
message("  DELETE_SIMPLICES  = ${DELETE_SIMPLICES}")
message("  DEBUG             = ${DEBUG}")
configure_file(
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h")
include_directories(${PROJECT_BINARY_DIR}) # because config.h file is in binary directory we need to include it


# --- sources and headers ---
set(SRC_LIST main.cpp tethex.cpp)
set(HDR_LIST tethex.h testing.h)


# --- Google Test ---
if(TESTING)
  message("ATTENTION! This is a version with testing features!")
  set(GTEST_DIR /home/artemiev/projects/gtest-1.6.0 CACHE PATH "Path to Google Test top level directory")
  find_library(GTEST_LIB gtest ${GTEST_DIR}/lib/.libs)
  if(GTEST_LIB)
    message("GTest library is ${GTEST_LIB}")
  else()
    message(FATAL_ERROR "GTest was not found!")
  endif()
  include_directories(${GTEST_DIR}/include)
endif()


# --- make ---
add_executable(${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
if(TESTING)
  target_link_libraries(${PROJECT_NAME} ${GTEST_LIB} pthread)
endif()


# --- packaging of the library ---
set(PACKAGE_FILE_NAME
    "${PROJECT_NAME}-${TETHEX_VERSION_MAJOR}.${TETHEX_VERSION_MINOR}.${TETHEX_VERSION_PATCH}")

add_custom_target(dist
                  COMMAND git archive --prefix=${PACKAGE_FILE_NAME}/ -o ${PROJECT_BINARY_DIR}/${PACKAGE_FILE_NAME}.tar.gz HEAD
                  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})


# --- autodocumentation using doxygen ---
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(
    "${PROJECT_SOURCE_DIR}/doxyfile.in"
    "${CMAKE_CURRENT_BINARY_DIR}/doxyfile")
  add_custom_target(doc COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxyfile)
else()
  message("WARNING: Doxygen is not found - Autodocumentation will not be created")
endif()

